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

int main()
{
	esc_data = 0;
	PIOS_Board_Init();

	PIOS_ADC_Config(1);
	
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

	// TODO: Move this into a PIOS_DELAY functi

	// TODO: Move this into a PIOS_DELAY function
	TIM_OCInitTypeDef tim_oc_init = {
		.TIM_OCMode = TIM_OCMode_PWM1,
		.TIM_OutputState = TIM_OutputState_Enable,
		.TIM_OutputNState = TIM_OutputNState_Disable,
		.TIM_Pulse = 0,
		.TIM_OCPolarity = TIM_OCPolarity_High,
		.TIM_OCNPolarity = TIM_OCPolarity_High,
		.TIM_OCIdleState = TIM_OCIdleState_Reset,
		.TIM_OCNIdleState = TIM_OCNIdleState_Reset,
	};
	TIM_OC1Init(TIM4, &tim_oc_init);
	TIM_ITConfig(TIM4, TIM_IT_CC1, ENABLE);  // Enabled by FSM

	NVIC_InitTypeDef NVIC_InitStructure;
	NVIC_InitStructure.NVIC_IRQChannel = TIM4_IRQn;
	NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGH;
	NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0;
	NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
	NVIC_Init(&NVIC_InitStructure);

	PIOS_LED_On(LED_GO);
	PIOS_LED_Off(LED_ERR);

	PIOS_ESC_Off();

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
		

		esc_process_static_fsm_rxn();

		// Serial interface: Process any incoming characters, and then process
		// any ongoing messages
		uint8_t c;
		if(PIOS_COM_ReceiveBuffer(PIOS_COM_DEBUG, &c, 1, 0) == 1)
			esc_serial_parse(c);
		esc_serial_process();
		
		if(esc_control.save_requested && esc_data->state == ESC_STATE_IDLE) {
			esc_control.save_requested = false;
			esc_settings_save(&config);
			// TODO: Send serial ack if succeeded or failed
		}
		
		if (esc_control.control_method == ESC_CONTROL_SERIAL) {
			esc_data->speed_setpoint = esc_control.serial_input;
		}
	}
	return 0;
}

#define MAX_RUNNING_FILTER 512
uint32_t DEMAG_BLANKING = 70;

#define MAX_CURRENT_FILTER 16
int32_t current_filter[MAX_CURRENT_FILTER];
int32_t current_filter_sum = 0;
int32_t current_filter_pointer = 0;

static int16_t diff_filter[MAX_RUNNING_FILTER];
static int32_t diff_filter_pointer = 0;
int32_t running_filter_length = 512;
static int32_t running_filter_sum = 0;

uint32_t calls_to_detect = 0;
uint32_t calls_to_last_detect = 0;

uint32_t samples_averaged;
int32_t threshold = -10;
int32_t divisor = 40;
#include "pios_adc_priv.h"
uint32_t detected;
uint32_t bad_flips;

uint32_t init_time;

// Index to look up filter lengths by RPM (incrementing by 128 rpm per entry
static const int32_t filter_length[] = {
	744, // ... 128
	372, // ... 256
	248, // ... 384
	186, // ... 512
	149, // ... 640
	124, // ... 768
	106, // ... 896
	93, // ... 1024
	83, // ... 1152
	64, // ... 1280
	58, // ... 1408
	52, // ... 1536
	47, // ... 1664
	43, // ... 1792
	40, // ... 1920
	37, // ... 2048
	34, // ... 2176
	31, // ... 2304
	29, // ... 2432
	27, // ... 2560
	25, // ... 2688
	24, // ... 2816
	22, // ... 2944
	21, // ... 3072
	20, // ... 3200
	19, // ... 3328
	18, // ... 3456
	18, // ... 3584
	17, // ... 3712
	17, // ... 3840
	17, // ... 3968
	16, // ... 4096
	16, // ... 4224
	16, // ... 4352
	15, // ... 4480
	15, // ... 4608
	15, // ... 4736
	14, // ... 4864
	14, // ... 4992
	14, // ... 5120
	14, // ... 5248
	14, // ... 5376
	13, // ... 5504
	13, // ... 5632
	13, // ... 5760
	13, // ... 5888
	13, // ... 6016
	12, // ... 6144
	12, // ... 6272
	12, // ... 6400
	12, // ... 6528
	12, // ... 6656
	12, // ... 6784
	12, // ... 6912
	11, // ... 7040
	11, // ... 7168
	11, // ... 7296
	11, // ... 7424
	11, // ... 7552
	10, // ... 7680
	10, // ... 7808
	10, // ... 7936
	9, // ... 8064
	9, // ... 8192
	9, // ... 8320
	9, // ... 8448
	9, // ... 8576
	9, // ... 8704
	9, // ... 8832
	8, // ... 8960
	8, // ... 9088
	8, // ... 9216
	8, // ... 9344
	7, // ... 9472
	7, // ... 9600
	7, // ... 9728
	7, // ... 9856
	6, // ... 9984
	6, // ... 10112
	6, // ... 10240
	6, // ... 10368
	6, // ... 10496
	6, // ... 10624
	6, // ... 10752
	5, // ... 10880
	5, // ... 11008
	5, // ... 11136
	5, // ... 11264
	5, // ... 11392
	5, // ... 11520
	5, // ... 11648
	4, // ... 11776
	4, // ... 11904
	4, // ... 12032
	4, // ... 12160
	4, // ... 12288
	4, // ... 12416
	3, // ... 12544
	3, // ... 12672
	3, // ... 12800
	3, // ... 12928
	2, // ... 13056
	2, // ... 13184
	2, // ... 13312
	2, // ... 13440
	2, // ... 13568
	2, // ... 13696
	2, // ... 13824
	2, // ... 13952
	2, // ... 14080
	2, // ... 14208
	2, // ... 14336
	2, // ... 14464
	2, // ... 14592
	2, // ... 14720
	2, // ... 14848
	2, // ... 14976
	2, // ... 15104
	2, // ... 15232
	2, // ... 15360
	2, // ... 15488
	2, // ... 15616
	2, // ... 15744
	2, // ... 15872
	2, // ... 16000
	2, // ... 16128
	2, // ... 16256
	2, // ... 16384
	2, // ... 16512
	2, // ... 16640
	2, // ... 16768
	2, // ... 16896
	2, // ... 17024
	2, // ... 17152
	2, // ... 17280
	2, // ... 17408
	2, // ... 17536
	2, // ... 17664
	2, // ... 17792
	2, // ... 17920
	1,
};	

uint8_t filter_length_scalar = 25;

const uint8_t undriven_pin[6] = {
	[ESC_STATE_AC] = 1,
	[ESC_STATE_CA] = 1,
	[ESC_STATE_AB] = 2,
	[ESC_STATE_BA] = 2,
	[ESC_STATE_BC] = 0,
	[ESC_STATE_CB] = 0
};

// This depends on the driving mode, assuming PWM high
const uint8_t high_pin[6] = {
	[ESC_STATE_AC] = 2,
	[ESC_STATE_CA] = 0,
	[ESC_STATE_AB] = 1,
	[ESC_STATE_BA] = 0,
	[ESC_STATE_BC] = 2,
	[ESC_STATE_CB] = 1
};

const bool pos[6] = {
	[ESC_STATE_AC] = true,
	[ESC_STATE_CA] = false,
	[ESC_STATE_AB] = false,
	[ESC_STATE_BA] = true,
	[ESC_STATE_BC] = false,
	[ESC_STATE_CB] = true
};

uint32_t adc_count = 0;
void DMA1_Channel1_IRQHandler(void)
{	
	adc_count++;	
	static enum pios_esc_state prev_state = ESC_STATE_AB;
//	static uint8_t overcurrent_count = 0;
	static int16_t * raw_buf;
	enum pios_esc_state curr_state;
	
	if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->full_flag /*DMA1_IT_TC1*/)) {	// whole double buffer filled
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[pios_adc_devs[0].dma_half_buffer_size];
		DMA_ClearFlag(pios_adc_devs[0].cfg->full_flag);
	}
	else if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->half_flag /*DMA1_IT_HT1*/)) {
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[0];
		DMA_ClearFlag(pios_adc_devs[0].cfg->half_flag);
	}
	else {
		// This should not happen, probably due to transfer errors
		DMA_ClearFlag(pios_adc_devs[0].cfg->dma.irq.flags /*DMA1_FLAG_GL1*/);
		return;
	}

	if(!esc_data)
		return;
	
	bad_flips += (raw_buf == pios_adc_devs[0].valid_data_buffer);
	raw_buf = (int16_t *) pios_adc_devs[0].valid_data_buffer;

	curr_state = PIOS_ESC_GetState();

	// If requested log data until buffer full
	if(esc_control.backbuffer_logging_status == ESC_LOGGING_CAPTURE) {
		// When T3 flag is off then it is sampled when driven
		uint16_t flags = 0;
		flags |= (TIM2->CR1 & TIM_CR1_DIR) ? 0x8000 : 0;
		flags |= (TIM3->CR1 & TIM_CR1_DIR) ? 0x4000 : 0;
		int16_t samples[4] = {raw_buf[0], raw_buf[1] | flags, raw_buf[2] | flags, raw_buf[3] | flags};
		if (esc_logger_put((uint8_t *) &samples[0], 8) != 0)
			esc_control.backbuffer_logging_status = ESC_LOGGING_FULL;
	}

	if( PIOS_DELAY_DiffuS(esc_data->last_swap_time) < DEMAG_BLANKING )
		return;
	
	// Smooth the estimate of current a bit 	
	int32_t current_measurement = raw_buf[0] * 8;// Convert to almost mA
	esc_data->current_ma += (current_measurement - esc_data->current_ma) / 50; // IIR filter
	

/*	if (esc_data->current > CURRENT_LIMIT && overcurrent_count > 4)
		esc_fsm_inject_event(ESC_EVENT_OVERCURRENT, 0);
	else if (esc_data->current > CURRENT_LIMIT)
			overcurrent_count ++;
	else overcurrent_count = 0; */

	   
	// If detected this commutation don't bother here
	if(curr_state == prev_state && esc_data->detected)
	   return;
	else if(curr_state != prev_state) {
		uint32_t timeval = PIOS_DELAY_GetRaw();
		
		prev_state = curr_state;
		calls_to_detect = 0;
		running_filter_sum = 0;
		diff_filter_pointer = 0;
		samples_averaged = 0;
		
		// Currently can get into a situation at high RPMs that causes filter to get
		// too short and then trigger lots of false ZCDs when it thinks its running at
		// 15k rpm
		uint32_t current_speed = esc_data->current_speed;
		if (current_speed > 8500)
			current_speed = 8500;
			
		if((current_speed >> 7) > NELEMENTS(filter_length))
		   running_filter_length = filter_length[NELEMENTS(filter_length)-1];
		else
		   running_filter_length = filter_length[esc_data->current_speed >> 7];

		running_filter_length = running_filter_length * filter_length_scalar / 100;
		
		if(running_filter_length >= MAX_RUNNING_FILTER)
			running_filter_length = MAX_RUNNING_FILTER;

		init_time = PIOS_DELAY_DiffuS(timeval);
	}

	int32_t last_swap_delay = PIOS_DELAY_DiffuS(esc_data->last_swap_time);
	if(last_swap_delay >= 0 && last_swap_delay < DEMAG_BLANKING)
		return;

	calls_to_detect++;
	
	int16_t undriven = raw_buf[1 + undriven_pin[curr_state]];
	int16_t high = raw_buf[1 + high_pin[curr_state]];
	int16_t ref = (raw_buf[1 + 0] + 
				   raw_buf[1 + 1] +
				   raw_buf[1 + 2]) / 3;
	int32_t diff;
	
	if(esc_data->duty_cycle > (PIOS_ESC_MAX_DUTYCYCLE * 0.05)) {
		// 127 / 27 is the 12.7 / 2.7 resistor divider
		// 3300 / 4096 is the mv / LSB
		uint32_t battery_mv = (high * PIOS_ESC_MAX_DUTYCYCLE * 127 * 3300) / (esc_data->duty_cycle * 27 * 4096);
		esc_data->battery_mv = (15 * esc_data->battery_mv + battery_mv) / 16;
	}
	
	if(pos[curr_state])
		diff = undriven - ref;
	else
		diff = ref - undriven;
	
	// Update running sum and history
	running_filter_sum += diff;
	// To avoid having to wipe the filter each commutation
	if(samples_averaged >= running_filter_length)
		running_filter_sum -= diff_filter[diff_filter_pointer];
	diff_filter[diff_filter_pointer] = diff;
	diff_filter_pointer++;
	if(diff_filter_pointer >= running_filter_length)
		diff_filter_pointer = 0;
	samples_averaged++;
	
	if(running_filter_sum > 0 && samples_averaged >= running_filter_length) {
		esc_data->detected = true;
		esc_fsm_inject_event(ESC_EVENT_ZCD, TIM4->CNT);
		calls_to_last_detect = calls_to_detect;
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
				PIOS_DELAY_WaitmS(1); //Count 1ms intervals in order to allow for possibility of watchdog
			}

			PIOS_LED_Off(LED_ERR);
			PIOS_LED_Off(PIOS_LED_HEARTBEAT);
			for(int i = 0 ; i < 250; i++) {
				PIOS_DELAY_WaitmS(1); //Count 1ms intervals in order to allow for possibility of watchdog
			}

		}
		PIOS_DELAY_WaitmS(1000);
	}
}

//TODO: Abstract out constants.  Need to know battery voltage too
//TODO: Other things to test for 
//      - impedance from motor(?)
//      - difference between high voltages
int32_t voltages[6][3];
void test_esc() {


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
	PIOS_DELAY_WaituS(3000);
	voltages[1][0] = PIOS_ADC_PinGet(1);
	voltages[1][1] = PIOS_ADC_PinGet(2);
	voltages[1][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_A_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(3000);
	voltages[0][0] = PIOS_ADC_PinGet(1);
	voltages[0][1] = PIOS_ADC_PinGet(2);
	voltages[0][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(3000);
	voltages[3][0] = PIOS_ADC_PinGet(1);
	voltages[3][1] = PIOS_ADC_PinGet(2);
	voltages[3][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(3000);
	voltages[2][0] = PIOS_ADC_PinGet(1);
	voltages[2][1] = PIOS_ADC_PinGet(2);
	voltages[2][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(3000);
	voltages[5][0] = PIOS_ADC_PinGet(1);
	voltages[5][1] = PIOS_ADC_PinGet(2);
	voltages[5][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(3000);
	voltages[4][0] = PIOS_ADC_PinGet(1);
	voltages[4][1] = PIOS_ADC_PinGet(2);
	voltages[4][2] = PIOS_ADC_PinGet(3);

	// If the particular phase isn't moving fet is dead
	if(voltages[0][0] < 1000) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_A_HIGH);
		panic(1);
	}
	if(voltages[1][0] > 30) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_A_LOW);
		panic(2);		
	}
	if(voltages[2][1] < 1000) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_B_HIGH);
		panic(3);
	}
	if(voltages[3][1] > 30){
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_B_LOW);
		panic(4);
	}	
	if(voltages[4][2] < 1000) {
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_C_HIGH);
		panic(5);
	}
	if(voltages[5][2] > 30){
		PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 10);
		PIOS_ESC_TestGate(ESC_C_LOW);
		panic(6);
	}
	// TODO: If other channels don't follow then motor lead bad
	
	PIOS_ESC_Off();

}

uint32_t bad_inputs;
void PIOS_TIM_4_irq_override();
extern void PIOS_DELAY_timeout();
void TIM4_IRQHandler(void) __attribute__ ((alias ("PIOS_TIM_4_irq_handler")));
static bool rising = false;

static uint16_t rise_value;
static uint16_t fall_value;
static uint16_t capture_value;

static void PIOS_TIM_4_irq_handler (void)
{
	static uint32_t last_input_update;
	
	if(TIM_GetITStatus(TIM4,TIM_IT_CC1)) {
		PIOS_DELAY_timeout();
		TIM_ClearITPendingBit(TIM4,TIM_IT_CC1);
	}
	
	if (TIM_GetITStatus(TIM4, TIM_IT_CC3)) {
		
		TIM_ClearITPendingBit(TIM4,TIM_IT_CC3);
		
		TIM_ICInitTypeDef TIM_ICInitStructure = {
			.TIM_ICPolarity = TIM_ICPolarity_Rising,
			.TIM_ICSelection = TIM_ICSelection_DirectTI,
			.TIM_ICPrescaler = TIM_ICPSC_DIV1,
			.TIM_ICFilter = 0x0,
		};
		
		if(rising) {
			rising = false;
			rise_value = TIM_GetCapture3(TIM4);
			
			/* Switch polarity of input capture */
			TIM_ICInitStructure.TIM_ICPolarity = TIM_ICPolarity_Falling;
			TIM_ICInitStructure.TIM_Channel = TIM_Channel_3;
			TIM_ICInit(TIM4, &TIM_ICInitStructure);
		} else {
			rising = true;
			fall_value = TIM_GetCapture3(TIM4);

			/* Switch polarity of input capture */
			TIM_ICInitStructure.TIM_ICPolarity = TIM_ICPolarity_Rising;
			TIM_ICInitStructure.TIM_Channel = TIM_Channel_3;
			TIM_ICInit(TIM4, &TIM_ICInitStructure);
			
			if (fall_value > rise_value) {
				capture_value = fall_value - rise_value;
			} else {
				capture_value = TIM4->ARR + fall_value - rise_value;
			}
			
		}
		
		// Dont process crazy values
		if(capture_value > config.PwmMin * 0.5 && capture_value < config.PwmMax * 1.1) {
			last_input_update = PIOS_DELAY_GetRaw();
			
			// Limit input range
			if(capture_value < config.PwmMin)
				capture_value = 0;
			else if (capture_value > config.PwmMax)
				capture_value = config.PwmMax;
			
			esc_control.pwm_input = capture_value;
			if(esc_control.control_method == ESC_CONTROL_PWM) {
				if (config.Mode == ESCSETTINGS_MODE_CLOSED) {
					uint32_t scaled_capture = config.RpmMin + (capture_value - config.PwmMin) * (config.RpmMax - config.RpmMin) / (config.PwmMax - config.PwmMin);
					esc_data->speed_setpoint = (capture_value < config.PwmMin) ? 0 : scaled_capture;
					esc_data->duty_cycle_setpoint = 0;
				} else if (config.Mode == ESCSETTINGS_MODE_OPEN){
					uint32_t scaled_capture = (capture_value - config.PwmMin) * PIOS_ESC_MAX_DUTYCYCLE / (config.PwmMax - config.PwmMin);
					esc_data->duty_cycle_setpoint = (capture_value < config.PwmMin) ? 0 : scaled_capture;
					esc_data->speed_setpoint = 0;
				} else {
					esc_data->duty_cycle_setpoint = 0;
					esc_data->speed_setpoint = 0;

				}
			}
		}
	} 
	
	if (TIM_GetITStatus(TIM4, TIM_IT_Update)) {
		if (PIOS_DELAY_DiffuS(last_input_update) > 100000) {
			esc_control.pwm_input = -1;
			if(esc_control.control_method == ESC_CONTROL_PWM) {
				esc_data->speed_setpoint = -1;
			}
		}
		TIM_ClearITPendingBit(TIM4,TIM_IT_Update);
	}
}

/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


