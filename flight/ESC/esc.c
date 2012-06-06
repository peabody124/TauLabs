/**
 ******************************************************************************
 * @addtogroup ESC esc
 * @brief The main ESC code
 *
 * @file       esc.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      INSGPS Test Program
 * @see        The GNU Public License (GPL) Version 3
 *
 *****************************************************************************/
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* OpenPilot Includes */
#include "pios.h"
#include "esc.h"

#include "fifo_buffer.h"
#include <pios_stm32.h>

#define CURRENT_LIMIT 4600

#define DOWNSAMPLING 1

//TODO: Check the ADC buffer pointer and make sure it isn't dropping swaps
//TODO: Check the time commutation is being scheduled, make sure it's the future
//TODO: Slave two timers together so in phase
//TODO: Ideally lock ADC and delay timers together to both
//TODO: Look into using TIM1
//TODO: Reenable watchdog and replace all PIOS_DELAY_WaitmS with something safe
//know the exact time of each sample and the PWM phase

//TODO: Measure battery voltage and normalize the feedforward model to be DC / Voltage

#define BACKBUFFER_ADC
//#define BACKBUFFER_ZCD
//#define BACKBUFFER_DIFF

/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

#define LED_GO  PIOS_LED_HEARTBEAT
#define LED_ERR PIOS_LED_ALARM

int16_t zero_current = 0;

const uint8_t dT = 1e6 / PIOS_ADC_RATE; // 6 uS per sample at 160k
float rate = 0;

static void test_esc();
static void panic(int diagnostic_code);
uint16_t pwm_duration ;
uint32_t counter = 0;

#define NUM_SETTLING_TIMES 20
uint32_t timer;
uint16_t timer_lower;
uint32_t step_period = 0x0080000;
uint32_t last_step = 0;
int16_t low_voltages[3];
int32_t avg_low_voltage;
struct esc_fsm_data * esc_data = 0;

/**
 * @brief ESC Main function
 */

uint32_t offs = 0;

// Major global variables
struct esc_control esc_control;
extern EscSettingsData config;
extern uint32_t pios_rcvr_group_map[1];
int main()
{
	esc_data = 0;
	PIOS_Board_Init();

	PIOS_ADC_Config(1);
	
	if (esc_settings_load(&config) != 0)
		esc_settings_defaults(&config);
	PIOS_ESC_SetPwmRate(config.PwmFreq);

	// TODO: Move this into an esc_control section
	esc_control.control_method = ESC_CONTROL_PWM;
	esc_control.serial_input = -1;
	esc_control.pwm_input = -1;
	esc_control.serial_logging_enabled = false;
	esc_control.save_requested = false;
	esc_control.backbuffer_logging_status = false;

	ADC_InitTypeDef ADC_InitStructure;
	ADC_StructInit(&ADC_InitStructure);
	ADC_InitStructure.ADC_Mode = ADC_Mode_RegSimult;
	ADC_InitStructure.ADC_ScanConvMode = ENABLE;
	ADC_InitStructure.ADC_ContinuousConvMode = DISABLE;
	ADC_InitStructure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T3_TRGO;
	ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
	ADC_InitStructure.ADC_NbrOfChannel = ((PIOS_ADC_NUM_CHANNELS + 1) >> 1);
	ADC_Init(ADC1, &ADC_InitStructure);
	ADC_Init(ADC2, &ADC_InitStructure);
	ADC_ExternalTrigConvCmd(ADC1, ENABLE);
	ADC_ExternalTrigConvCmd(ADC2, ENABLE);

	PIOS_LED_On(LED_GO);
	PIOS_LED_Off(LED_ERR);

	PIOS_ESC_Off();
	PIOS_ESC_SetDirection(config.Direction == ESCSETTINGS_DIRECTION_FORWARD ? ESC_FORWARD : ESC_BACKWARD);

	esc_serial_init();
	
	test_esc();
	
	// Blink LED briefly once passed	
	PIOS_LED_Off(0);
	PIOS_LED_Off(1);
	PIOS_DELAY_WaitmS(250);
	PIOS_LED_On(0);
	PIOS_LED_On(1);
	PIOS_DELAY_WaitmS(500);
	PIOS_LED_Off(0);
	PIOS_LED_Off(1);
	PIOS_DELAY_WaitmS(250);
	
	esc_data = esc_fsm_init();
	esc_data->speed_setpoint = -1;
	esc_data->duty_cycle_setpoint = -1;

	PIOS_ADC_StartDma();
	
	counter = 0;
	uint32_t timeval = PIOS_DELAY_GetRaw();
	uint32_t ms_count = 0;
	while(1) {
		counter++;
		
		if(PIOS_DELAY_DiffuS(timeval) > 1000) {
			ms_count++;
			timeval = PIOS_DELAY_GetRaw();
			// Flash LED every 1024 ms
			if((ms_count & 0x000007ff) == 0x400) {
				PIOS_LED_Toggle(0);
			}

			if (esc_control.serial_logging_enabled) {
				uint16_t send_buffer[8] = {0xff00, (ms_count & 0x0000ffff), (ms_count & 0xffff0000) >> 16, esc_data->current_speed, esc_data->speed_setpoint, esc_data->duty_cycle};
				PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, (uint8_t *) send_buffer, sizeof(send_buffer));
			}
		}
		
		// 1 based indexing on channels
		uint16_t capture_value = PIOS_RCVR_Read(pios_rcvr_group_map[0],1);
		esc_control.pwm_input = capture_value;
		if(esc_control.control_method == ESC_CONTROL_PWM) {
			if (capture_value == (uint16_t) PIOS_RCVR_INVALID || capture_value == (uint16_t) PIOS_RCVR_NODRIVER || (uint16_t)	capture_value == PIOS_RCVR_TIMEOUT) {
				esc_data->speed_setpoint = -1;
				esc_data->duty_cycle_setpoint = -1;
			}
			else if (config.Mode == ESCSETTINGS_MODE_CLOSED) {
				uint32_t scaled_capture = config.RpmMin + (capture_value - config.PwmMin) * (config.RpmMax - config.RpmMin) / (config.PwmMax - config.PwmMin);
				esc_data->speed_setpoint = (capture_value < config.PwmMin) ? 0 : scaled_capture;
				esc_data->duty_cycle_setpoint = -1;
			} else if (config.Mode == ESCSETTINGS_MODE_OPEN){
				uint32_t scaled_capture = (capture_value - config.PwmMin) * PIOS_ESC_MAX_DUTYCYCLE / (config.PwmMax - config.PwmMin);
				esc_data->duty_cycle_setpoint = (capture_value < config.PwmMin) ? 0 : scaled_capture;
				esc_data->speed_setpoint = -1;
			} else { // Invalid mode specified
				esc_data->speed_setpoint = -1;
				esc_data->duty_cycle_setpoint = -1;
			}
		}

		esc_process_static_fsm_rxn();

		// Serial interface: Process any incoming characters, and then process
		// any ongoing messages
		uint8_t c;
		if(PIOS_COM_ReceiveBuffer(PIOS_COM_DEBUG, &c, 1, 0) == 1)
			esc_serial_parse(c);
		esc_serial_process();
		
		if(esc_control.save_requested && (esc_data->state == ESC_STATE_IDLE || esc_data->state == ESC_STATE_WAIT_FOR_ARM)) {
			esc_control.save_requested = false;
			int8_t ret_val = esc_settings_save(&config);
			PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, (uint8_t *) &ret_val, sizeof(ret_val));
			// TODO: Send serial ack if succeeded or failed
		}
		
		if (esc_control.control_method == ESC_CONTROL_SERIAL) {
			if (config.Mode == ESCSETTINGS_MODE_CLOSED)
				esc_data->speed_setpoint = esc_control.serial_input;
			else if (config.Mode == ESCSETTINGS_MODE_OPEN)
				esc_data->duty_cycle_setpoint = PIOS_ESC_MAX_DUTYCYCLE * esc_control.serial_input / 10000;
		}
	}
	return 0;
}

void panic_wait(uint32_t delay_ms)
{
	// Serial interface: Process any incoming characters, and then process
	// any ongoing messages
	uint32_t in_time = PIOS_DELAY_GetRaw();
	
	while(PIOS_DELAY_DiffuS(in_time) < (delay_ms * 1000)) {
		uint8_t c;
		if(PIOS_COM_ReceiveBuffer(PIOS_COM_DEBUG, &c, 1, 0) == 1)
			esc_serial_parse(c);
		esc_serial_process();
	}
}

/* INS functions */
void panic(int diagnostic_code)
{
	// Turn off error LED.
	PIOS_LED_Off(LED_ERR);
	while(1) {
		for(int i=0; i<diagnostic_code; i++)
		{
			PIOS_LED_On(LED_ERR);
			PIOS_LED_On(PIOS_LED_HEARTBEAT);
			for(int i = 0 ; i < 250; i++) {
				panic_wait(1); //Count 1ms intervals in order to allow for possibility of watchdog
			}

			PIOS_LED_Off(LED_ERR);
			PIOS_LED_Off(PIOS_LED_HEARTBEAT);
			for(int i = 0 ; i < 250; i++) {
				panic_wait(1); //Count 1ms intervals in order to allow for possibility of watchdog
			}

		}
		panic_wait(1000);
	}
}

//TODO: Abstract out constants.  Need to know battery voltage too
//TODO: Other things to test for 
//      - impedance from motor(?)
//      - difference between high voltages
int16_t voltages[6][3];
static void get_voltages()
{
	PIOS_ESC_Off();
	for(int i = 0; i < 150; i++) {
		PIOS_DELAY_WaitmS(1);
	}
	
	zero_current = PIOS_ADC_PinGet(0);
	
	PIOS_ESC_Arm();
	
	PIOS_ESC_TestGate(ESC_A_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[0] = PIOS_ADC_PinGet(1);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[1] = PIOS_ADC_PinGet(2);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[2] = PIOS_ADC_PinGet(3);
	avg_low_voltage = low_voltages[0] + low_voltages[1] + low_voltages[2];
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_A_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[1][0] = PIOS_ADC_PinGet(1);
	voltages[1][1] = PIOS_ADC_PinGet(2);
	voltages[1][2] = PIOS_ADC_PinGet(3);
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_A_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[0][0] = PIOS_ADC_PinGet(1);
	voltages[0][1] = PIOS_ADC_PinGet(2);
	voltages[0][2] = PIOS_ADC_PinGet(3);
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[3][0] = PIOS_ADC_PinGet(1);
	voltages[3][1] = PIOS_ADC_PinGet(2);
	voltages[3][2] = PIOS_ADC_PinGet(3);
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[2][0] = PIOS_ADC_PinGet(1);
	voltages[2][1] = PIOS_ADC_PinGet(2);
	voltages[2][2] = PIOS_ADC_PinGet(3);
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[5][0] = PIOS_ADC_PinGet(1);
	voltages[5][1] = PIOS_ADC_PinGet(2);
	voltages[5][2] = PIOS_ADC_PinGet(3);
	
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(250);
	voltages[4][0] = PIOS_ADC_PinGet(1);
	voltages[4][1] = PIOS_ADC_PinGet(2);
	voltages[4][2] = PIOS_ADC_PinGet(3);
}

const float scale = 3.3 / 4096 * 12.7 / 2.7;

float max_lead_difference(int16_t voltages[3])
{
	float d1, d2, d3;
	d1 = scale * fabs(voltages[0] - voltages[1]);
	d2 = scale * fabs(voltages[0] - voltages[2]);
	d3 = scale * fabs(voltages[1] - voltages[2]);
	
	float max = (d1 > d2) ? d1 : d2;
	max = (max > d3) ? max : d3;
	
	return max;
}

/**
 * Analyze the voltages from the initial power up test to determine
 * if a motor is present.  This should be improved to look at current
 * dynamics and extract the R / L values for the motor
 */
static bool check_motor()
{
	const float max_diff = 0.3;
	
	bool present = (max_lead_difference(voltages[0]) < max_diff) &&
		(max_lead_difference(voltages[1]) < max_diff) &&
		(max_lead_difference(voltages[2]) < max_diff) &&
		(max_lead_difference(voltages[3]) < max_diff) &&
		(max_lead_difference(voltages[4]) < max_diff) &&
		(max_lead_difference(voltages[5]) < max_diff) &&
		(max_lead_difference(voltages[6]) < max_diff);

	// For now just indicate this with long LED pulse
	if (!present) {
		PIOS_LED_On(LED_ERR);
		PIOS_DELAY_WaitmS(1000);
		PIOS_LED_Off(LED_ERR);
	}

	return present;
}

void test_esc() 
{
	const float max_low = 0.1;
	const bool testing = false;
	get_voltages();

	float battery = (voltages[0][0] + voltages[2][1] + voltages[4][2]) / 3 * scale;
	if (battery < 7) // Don't run if battery is too low
		panic(8);
	
	// Allow up to 5% difference on the channel values
	const float max_diff = 0.05 * battery;

	// If the particular phase isn't moving fet is dead
	if(fabs(voltages[0][0] * scale - battery) > max_diff) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_A_HIGH);
		panic(1);
	}
	if(voltages[1][0] * scale > max_low) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_A_LOW);
		panic(2);		
	}
	if(fabs(voltages[2][1] * scale - battery) > max_diff) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_B_HIGH);
		panic(3);
	}
	if(voltages[3][1] * scale > max_low){
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_B_LOW);
		panic(4);
	}	
	if(fabs(voltages[4][2] * scale - battery) > max_diff) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_C_HIGH);
		panic(5);
	}
	if(voltages[5][2] * scale > max_low){
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		if (testing) PIOS_ESC_TestGate(ESC_C_LOW);
		panic(6);
	}
	
	PIOS_ESC_Off();

	check_motor();
}
