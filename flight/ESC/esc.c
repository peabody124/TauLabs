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
#include "esc_fsm.h"
#include "fifo_buffer.h"
#include <pios_stm32.h>

//TODO: Check the ADC buffer pointer and make sure it isn't dropping swaps
//TODO: Check the time commutation is being scheduled, make sure it's the future
//TODO: Slave two timers together so in phase
//TODO: Ideally lock ADC and delay timers together to both
//TODO: Look into using TIM1
//know the exact time of each sample and the PWM phase

#define BACKBUFFER
#define BACKBUFFER_ZCD

/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

void adc_callback(float * buffer);

#define DOWNSAMPLING 6

#ifdef BACKBUFFER
uint16_t back_buf[8096];
uint16_t back_buf_point = 0;
#endif

#define LED_ERR LED1
#define LED_GO  LED2
#define LED_MSG LED3

// Tuning settings for control
float state_offset[6] = {40, 40, 40, 40, 40, 40};

uint16_t zero_current = 0;

volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin;
volatile bool pos;

const uint8_t dT = 1e6 / PIOS_ADC_RATE; // 6 uS per sample at 160k
float rate = 0;

static void test_esc();
static void panic(int diagnostic_code);
uint16_t pwm_duration ;
uint32_t counter = 0;

#define NUM_SETTLING_TIMES 20
uint32_t timer;
uint16_t timer_lower;
uint32_t step_period = 0x0100000;
uint32_t last_step = 0;
uint32_t settling_time[NUM_SETTLING_TIMES];
uint8_t settling_time_pointer = 0;
uint8_t settled = 1;

struct esc_fsm_data * esc_data = 0;

/**
 * @brief ESC Main function
 */
int main()
{
	esc_data = 0;
	PIOS_Board_Init();

	PIOS_ADC_Config(DOWNSAMPLING);
	PIOS_ADC_SetCallback(adc_callback);

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
	TIM_ITConfig(TIM4, TIM_IT_CC1, DISABLE);  // Enabled by FSM

	NVIC_InitTypeDef NVIC_InitStructure;
	NVIC_InitStructure.NVIC_IRQChannel = TIM4_IRQn;
	NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGH;
	NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0;
	NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
	NVIC_Init(&NVIC_InitStructure);

	// This pull up all the ADC voltages so the BEMF when at -0.7V
	// is still positive
	PIOS_GPIO_Enable(0);
	PIOS_GPIO_Off(0);

	PIOS_LED_Off(LED_ERR);
	PIOS_LED_On(LED_GO);
	PIOS_LED_On(LED_MSG);

	test_esc();
	esc_data = esc_fsm_init();
	esc_data->speed_setpoint = 2500;

	while(1) {
		counter++;

		// Create large size timer
		timer_lower = PIOS_DELAY_GetuS();
		if(timer_lower < (timer & 0x0000ffff))
			timer += 0x00010000;
		timer = (timer & 0xffff0000) | timer_lower;

		if((timer - last_step) > step_period) {
			last_step = timer;
			esc_data->speed_setpoint = (esc_data->speed_setpoint == 4500) ? 3500 : 4500;
			settled = 0;
		}

		if(!settled && abs(esc_data->speed_setpoint - esc_data->current_speed) < 100) {
			settled = 1;
			settling_time[settling_time_pointer++] = timer - last_step;
			if(settling_time_pointer > NUM_SETTLING_TIMES)
				settling_time_pointer = 0;
		}

		esc_process_static_fsm_rxn();

		switch(PIOS_ESC_GetState()) {
			case ESC_STATE_AC:
				low_pin = 0;
				high_pin = 2;
				undriven_pin = 1;
				pos = true;
				break;
			case ESC_STATE_CA:
				low_pin = 2;
				high_pin = 0;
				undriven_pin = 1;
				pos = false;
				break;
			case ESC_STATE_AB:
				low_pin = 0;
				high_pin = 1;
				undriven_pin = 2;
				pos = false;
				break;
			case ESC_STATE_BA:
				low_pin = 1;
				high_pin = 0;
				undriven_pin = 2;
				pos = true;
				break;
			case ESC_STATE_BC:
				low_pin = 1;
				high_pin = 2;
				undriven_pin = 0;
				pos = false;
				break;
			case ESC_STATE_CB:
				low_pin = 2;
				high_pin = 1;
				undriven_pin = 0;
				pos = true;
				break;
			default:
				PIOS_ESC_Off();
		}
	}
	return 0;
}

// When driving both legs to ground mid point is 580
#define ZERO_POINT 290
#define CROSS_ALPHA 0.5
#define MIN_PRE_COUNT  2
#define MIN_POST_COUNT 1
#define DEMAG_BLANKING 60
uint16_t exceed_count = 0;

uint32_t calls = 0;
uint32_t detected = 0;

void adc_callback(float * buffer)
{
	calls ++;

	int16_t * raw_buf = PIOS_ADC_GetRawBuffer();
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	enum pios_esc_state curr_state;
	static uint16_t below_time;
	static uint16_t pre_count = 0;
	static uint16_t post_count = 0;
	static int16_t running_avg;
	bool first;

	curr_state = PIOS_ESC_GetState();

#ifdef BACKBUFFER_ADC
	// Debugging code - keep a buffer of old ADC values
	if((back_buf_point + 3 * DOWNSAMPLING + 3) > (sizeof(back_buf) / sizeof(back_buf[0])))
		back_buf_point = 0;
	back_buf[back_buf_point++] = 0xFFFF;
	back_buf[back_buf_point++] = curr_state;
	back_buf[back_buf_point++] = 0;
	back_buf[back_buf_point++] = 0;
	for(int i = 0; i < DOWNSAMPLING; i++) {
		back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1];
		back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 2];
		back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 3];
	}
#endif

	// Wait for blanking
	if(esc_data && PIOS_DELAY_DiffuS(esc_data->last_swap_time) < DEMAG_BLANKING)
		return;

	if(esc_data) {
		// Smooth the estimate of current a bit
		int32_t this_current = 0;
		for(int i = 0; i < DOWNSAMPLING; i++)
			this_current += raw_buf[PIOS_ADC_NUM_CHANNELS * i + 0];
		this_current /= DOWNSAMPLING;
		this_current -= zero_current;
		if(this_current < 0)
			this_current = 0;
		esc_data->current += (this_current - esc_data->current) * 0.01;
	}

	uint16_t enter_time = PIOS_DELAY_GetuS();

	// Already detected for this stage
	if(esc_data && esc_data->detected)
		return;

	// If detected this commutation don't bother here
	if(curr_state != prev_state) {
		prev_state = curr_state;
		pre_count = 0;
		post_count = 0;
		first = true;
	} else {
		first = false;
	}


	switch(PIOS_ESC_GetMode()) {
		case ESC_MODE_LOW_ON_PWM_HIGH:
			// Doesn't work quite right yet
			for(int i = 0; i < DOWNSAMPLING; i++) {
				int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
				int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
				int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
				int16_t ref = (high + low) / 2;
				int16_t diff;

				if(pos)
					diff = undriven - ref - ZERO_POINT;
				else
					diff = ref - undriven - ZERO_POINT;

				diff *= 10; // scale up to increase fixed point precision

				if(first)
					running_avg = diff;
				else
					running_avg = (running_avg + diff) / 2; //(1-CROSS_ALPHA) * running_avg + CROSS_ALPHA * diff;

				diff = running_avg;

				if(diff < 0) {
					pre_count++;
					// Keep setting so we store the time of zero crossing
					below_time = enter_time - dT * (DOWNSAMPLING - i);
				}
				if(diff > 0 && pre_count > MIN_PRE_COUNT) {
					post_count++;
				}
				if(diff > 0 && post_count >= MIN_POST_COUNT) {

					detected ++;
					esc_data->detected = true;
					esc_fsm_inject_event(ESC_EVENT_ZCD, below_time);
#ifdef BACKBUFFER_ZCD
					back_buf[back_buf_point++] = below_time;
					back_buf[back_buf_point++] = esc_data->speed_setpoint;
					if(back_buf_point > (sizeof(back_buf) / sizeof(back_buf[0])))
						back_buf_point = 0;
#endif /* BACKBUFFER_ZCD */
					break;
				}
			}
			break;
		default:
			PIOS_ESC_Off();
			break;
	}
}

/* INS functions */
void panic(int diagnostic_code)
{
	PIOS_ESC_Off();
	// Polarity backwards
	PIOS_LED_On(LED_ERR);
	while(1) {
		for(int i=0; i<diagnostic_code; i++)
		{
			PIOS_LED_Toggle(LED_ERR);
			PIOS_DELAY_WaitmS(250);
			PIOS_LED_Toggle(LED_ERR);
			PIOS_DELAY_WaitmS(250);
		}
		PIOS_DELAY_WaitmS(1000);
	}
}


//TODO: Abstract out constants.  Need to know battery voltage too
void test_esc() {
	int32_t voltages[6][3];

	PIOS_DELAY_WaitmS(150);

	zero_current = PIOS_ADC_PinGet(0);

	PIOS_ESC_Arm();

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_A_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[1][0] = PIOS_ADC_PinGet(1);
	voltages[1][1] = PIOS_ADC_PinGet(2);
	voltages[1][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_A_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[0][0] = PIOS_ADC_PinGet(1);
	voltages[0][1] = PIOS_ADC_PinGet(2);
	voltages[0][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[3][0] = PIOS_ADC_PinGet(1);
	voltages[3][1] = PIOS_ADC_PinGet(2);
	voltages[3][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_B_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[2][0] = PIOS_ADC_PinGet(1);
	voltages[2][1] = PIOS_ADC_PinGet(2);
	voltages[2][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[5][0] = PIOS_ADC_PinGet(1);
	voltages[5][1] = PIOS_ADC_PinGet(2);
	voltages[5][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(0.5);
	PIOS_ESC_TestGate(ESC_C_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(1);
	PIOS_DELAY_WaituS(100);
	voltages[4][0] = PIOS_ADC_PinGet(1);
	voltages[4][1] = PIOS_ADC_PinGet(2);
	voltages[4][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_Off();
	// If the particular phase isn't moving fet is dead
	if(voltages[0][0] < 1000)
		panic(1);
	if(voltages[1][0] > 700)
		panic(2);
	if(voltages[2][1] < 1000)
		panic(2);
	if(voltages[3][1] > 700)
		panic(3);
	if(voltages[4][2] < 1000)
		panic(4);
	if(voltages[5][2] > 700)
		panic(5);

	// TODO: If other channels don't follow then motor lead bad
}

/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


