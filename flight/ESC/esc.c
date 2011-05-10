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
#include "fifo_buffer.h"
#include <pios_stm32.h>


//TODO: Check the ADC buffer pointer and make sure it isn't dropping swaps
//TODO: Check the time commutation is being scheduled, make sure it's the future
//TODO: Slave two timers together so in phase
//TODO: Ideally lock ADC and delay timers together to both 
//TODO: Look into using TIM1 
//know the exact time of each sample and the PWM phase

/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

void adc_callback(float * buffer);

#define DOWNSAMPLING 10

uint16_t back_buf[8096];
uint16_t back_buf_point = 0;

// A message from the ADC to say a zero crossing was detected
struct zerocrossing_message {
	enum pios_esc_state state;
	volatile uint16_t time;
	volatile bool read;
} message;

struct zerocrossing_stats {
	uint32_t consecutive_skipped;
	uint32_t consecutive_detected;
	uint16_t interval;
	uint16_t smoothed_interval;
	uint16_t latency[6];
} zerocrossing_stats;

bool closed_loop_updated = false;

// Tuning settings for control
float state_offset[6] = {40, 40, 40, 40, 40, 40};
float commutation_phase = 0.45;
float kp = 0.0001;
float ki = 0.0000001;
float kff = 5e-5;
float kff2 = 0.06;
float accum = 0.0;
float ilim = 0.5;

float max_dc_change = 0.001;

#define RPM_TO_US(x) (1e6 * 60 / (x * COMMUTATIONS_PER_ROT) )
#define US_TO_RPM(x) RPM_TO_US(x)
#define COMMUTATIONS_PER_ROT (7*6)
#define MIN_DC 0.01
#define MAX_DC 0.7
int16_t initial_startup_speed = 150;
int16_t final_startup_speed = 1000;
int16_t current_speed;
bool closed_loop = false;
int16_t desired_rpm = 2000;
volatile uint16_t swap_time;
volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin; 
volatile bool pos;
volatile float dc = 0.18;

volatile bool commutation_detected = false;
volatile bool commutated = false;
volatile uint16_t checks = 0;
uint16_t consecutive_nondetects = 0;

const uint8_t dT = 1e6 / PIOS_ADC_RATE; // 6 uS per sample at 160k
float rate = 0;

uint16_t swap_diff;

uint32_t idle_count = 0;
uint32_t last_idle_count;
bool schedule_next = false;

void process_message(struct zerocrossing_message * msg);
void commutate();

enum init_state {
	INIT_GRAB = 0,
	INIT_ACCEL = 1,
	INIT_WAIT = 2,
	INIT_FAIL = 3
} init_state;

void PIOS_DELAY_timeout() {
	TIM_ClearITPendingBit(TIM4,TIM_IT_CC1);
	TIM_ClearFlag(TIM4,TIM_IT_CC1);
	commutate();
}

void schedule_commutation(uint16_t time) 
{
	TIM_SetCompare1(TIM4, time);
}

/**
 * @brief ESC Main function
 */
int main()
{
	current_speed = initial_startup_speed;
	message.read = true;
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
	TIM_ITConfig(TIM4, TIM_IT_CC1, ENABLE);
	
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
		
	PIOS_LED_Off(LED1);
	PIOS_LED_On(LED2);
	PIOS_LED_On(LED3);
	
	PIOS_ESC_SetDutyCycle(dc);
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);
	PIOS_ESC_Arm();
	
	init_state = INIT_GRAB;
	closed_loop = false;
	
	while(1) {
		idle_count++;
		if(commutated) {
			last_idle_count = idle_count;
			idle_count = 0;
			commutated = false;
			
			if(closed_loop) {	
				// Turn err light off
				PIOS_LED_On(LED2);
				
				if(!commutation_detected) {
					consecutive_nondetects++;
					PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "- %u\n", checks);
				} else 
					consecutive_nondetects = 0;
				
				if(!closed_loop_updated) {
					PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "*");
					static uint16_t fail_count = 0;
					if(fail_count++ > 0) {
						PIOS_ESC_Off();
						TIM_ITConfig(TIM4, TIM_IT_CC1, DISABLE);
					}					
				}
				closed_loop_updated = false;
				
				if(consecutive_nondetects > 50) {
					PIOS_ESC_Off();
					TIM_ITConfig(TIM4, TIM_IT_CC1, DISABLE);
				}
				
				// This is a fall back.  Should get rescheduled by zero crossing detection.
				if(schedule_next)
					schedule_commutation(swap_time + zerocrossing_stats.smoothed_interval);
				else
					schedule_commutation(swap_time + 7 * zerocrossing_stats.smoothed_interval);
												
				int16_t error = desired_rpm - US_TO_RPM(zerocrossing_stats.smoothed_interval);
				
				accum += ki * error;
				if(accum > ilim)
					accum = ilim;
				if(accum < -ilim)
					accum = -ilim;
				
				float new_dc = desired_rpm * kff - kff2 + kp * error + accum;
				
				if((new_dc - dc) > max_dc_change)
					dc += max_dc_change;
				else if((new_dc - dc) < -max_dc_change)
					dc -= max_dc_change;
				else
					dc = new_dc;
				
				if(dc < MIN_DC)
					dc = MIN_DC;
				if(dc > MAX_DC)
					dc = MAX_DC;
				PIOS_ESC_SetDutyCycle(dc); 

				
				/*uint8_t state = PIOS_ESC_GetState();
				uint16_t desired_latency = zerocrossing_stats.smoothed_interval * (1 - commutation_phase);
				if(zerocrossing_stats.latency[state] > desired_latency)
					state_offset[state] -= 0.01;
				else if (zerocrossing_stats.latency[state] < desired_latency)
					state_offset[state] += 0.01;    */

				/*static uint16_t count = 0;
				if((count++ % 100) == 0)
					PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG,"%u\n", zerocrossing_stats.smoothed_interval); */
				
				//PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u\n", swap_diff);
				
				
			} else {
				// Turn err light off
				PIOS_LED_Off(LED2);
				
				static uint16_t init_counter;
				uint16_t delay = 0;
				
				// Simple startup state machine.  Needs constants removing
				switch(init_state) {
					case INIT_GRAB:
						PIOS_ESC_SetState(0);
						current_speed = initial_startup_speed;
						delay = 30000; // hold in that position for 10 ms
						PIOS_ESC_SetDutyCycle(0.2); // TODO: Current limit
						dc = 0.12;
						init_state = INIT_ACCEL;
						break;
					case INIT_ACCEL:
						if(current_speed < final_startup_speed) 
							current_speed+=2.0;
						else
							init_state = INIT_WAIT;
						init_counter = 0;
						delay = RPM_TO_US(current_speed);
						PIOS_ESC_SetDutyCycle(dc);
						break;
					case INIT_WAIT:
						dc = 0.1 + (dc - 0.1) * 0.999;
//						accum = dc;
						PIOS_ESC_SetDutyCycle(dc);
						delay = RPM_TO_US(current_speed);
						if(init_counter++ > 2000)
							init_state = INIT_FAIL;
						break;
					case INIT_FAIL:
						PIOS_ESC_Off();
						break;
				}
				
				schedule_commutation(PIOS_DELAY_GetuS() + delay);
				
			}			
			
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
			commutation_detected = false;
			checks = 0;
		}			
	}
	return 0;
}	

void commutate() 
{
	//PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", (next - swap_time), (uint32_t) (dc * 10000));
	uint16_t this_swap_time = PIOS_DELAY_GetuS();
	swap_diff = (this_swap_time - swap_time);
	swap_time = this_swap_time;
	
	PIOS_ESC_NextState();
	
	commutated = true;
}

uint16_t enter_time = 0;

/**
 * @brief Process any message from the zero crossing detection
 * @param[in] msg The message containing the state and time it occurred in
 */
void process_message(struct zerocrossing_message * msg) 
{
	bool skipped = false;
	static bool prev_skipped = false;
	static uint16_t last_time;
	static enum pios_esc_state prev_state;

	// Don't reprocess any read messages
	if(msg->read)
		return;	
	msg->read = true;
	
	// Sanity check
	if(msg->state != PIOS_ESC_GetState()) {
		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "."); //---------->Oh shit\n");
		return;
	}

	// Check for any skipped states
	switch(msg->state) {
		case 1:
			if(prev_state != 3)
				skipped = true;
			break;
		case 0:
			if(prev_state != 1)
				skipped = true;
			break;
		case 5:
			if(prev_state != 0)
				skipped = true;
			break;
		case 4:
			if(prev_state != 5)
				skipped = true;
			break;
		case 2:
			if(prev_state != 4)
				skipped = true;
			break;
		case 3:
			if(prev_state != 2)
				skipped = true;
			break;
		default:
			break;
	}
	
	if(!skipped) {
		PIOS_LED_Off(LED3);

		zerocrossing_stats.consecutive_skipped = 0;
		zerocrossing_stats.consecutive_detected++;
	} else {
		PIOS_LED_On(LED3);

		zerocrossing_stats.consecutive_skipped++;
		zerocrossing_stats.consecutive_detected = 0;
	}		
		
	// If meant to be working and missed a bunch
	if(closed_loop && zerocrossing_stats.consecutive_skipped > 50)
		PIOS_ESC_Off();

	// Compute interval since last zcd
	zerocrossing_stats.interval = msg->time - last_time;
	zerocrossing_stats.latency[msg->state] = zerocrossing_stats.latency[msg->state] * 0.9 + (msg->time - swap_time) * 0.1;
	last_time = message.time;
	
	if(skipped && closed_loop)
		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", msg->state, zerocrossing_stats.smoothed_interval);

	// If decent interval use it to update estimate of speed
	if(!skipped && !prev_skipped && (zerocrossing_stats.interval < 10000))
		zerocrossing_stats.smoothed_interval = 0.95 * zerocrossing_stats.smoothed_interval + 0.05 * zerocrossing_stats.interval;

	if(zerocrossing_stats.consecutive_detected > 82) 
		closed_loop = true;
							
	if(closed_loop) {
		// TOOD: This logic shouldn't stay here
		closed_loop_updated = true;
		if(PIOS_DELAY_DiffuS(msg->time + zerocrossing_stats.smoothed_interval * commutation_phase) >= -2) {
			// Fallback
			PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "# %u %u %u %u %u %u %u\n", msg->time, swap_time, PIOS_DELAY_GetuS(), enter_time, zerocrossing_stats.smoothed_interval, msg->state, commutated);
			schedule_commutation(PIOS_DELAY_GetuS() + 2);
		} else {
			// Keep all math uint16_t to deal with timer wrap around
			uint16_t zcd_time = msg->time + zerocrossing_stats.smoothed_interval * commutation_phase;
			uint16_t zcd_latency = zcd_time - swap_time;
			uint16_t predicted_latency = zerocrossing_stats.smoothed_interval;
			uint16_t schedule_time = swap_time + (predicted_latency*0 + 2*zcd_latency) / 2;
	
			schedule_commutation(schedule_time);
		}
	}

	
	prev_skipped = skipped;
	prev_state = msg->state;
}

// When driving both legs to ground mid point is 580
#define MID_POINT 580
#define DIODE_LOW 460
#define MIN_PRE_COUNT  1
#define MIN_POST_COUNT 1
#define DEMAG_BLANKING 20
int16_t big_val;
uint16_t exceed_count = 0;

volatile bool adc_irq_status = false;
volatile uint32_t adc_irq_collisions = 0;

void adc_callback(float * buffer) 
{
	if(adc_irq_status)
		adc_irq_collisions++;
	adc_irq_status = true;
		
	int16_t * raw_buf = PIOS_ADC_GetRawBuffer();
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	enum pios_esc_state curr_state;	
	static uint16_t below_time;
	static uint16_t pre_count = 0;
	static uint16_t post_count = 0;
	static int16_t running_avg;

	// Commutation detection, assuming mode is ESC_MODE_LOW_ON_PWM_BOTH
	/* uint16_t */ enter_time = PIOS_DELAY_GetuS();
	// Process ADC from undriven leg
	// TODO: Make these variables be updated when a state transition occurs

	// Wait for blanking
	if(PIOS_DELAY_DiffuS(swap_time)  < DEMAG_BLANKING)
		goto adc_done;
	
	if(commutation_detected)
		goto adc_done;
	
	// If detected this commutation don't bother here
	// TODO:  disable IRQ for efficiency
	curr_state = PIOS_ESC_GetState();

	if(curr_state != prev_state) {
		prev_state = curr_state;

		pre_count = 0;
		post_count = 0;
		running_avg = 0;
	}

	checks++;
	switch(PIOS_ESC_GetMode()) {
		case ESC_MODE_LOW_ON_PWM_HIGH:
			// Doesn't work quite right yet
			if((back_buf_point + 3 * DOWNSAMPLING + 2) > (sizeof(back_buf) / 2))
				back_buf_point = 0;
			back_buf[back_buf_point++] = 0xFFFF;
			back_buf[back_buf_point++] = curr_state;
			for(int i = 0; i < DOWNSAMPLING; i++) {
				int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
				int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
				int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
				int16_t ref = (high + MID_POINT) / 2; 
				int16_t diff;
				
				back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1];
				back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 2];
				back_buf[back_buf_point++] = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 3];
								
				if(pos) {						
					diff = undriven - MID_POINT - state_offset[curr_state]; 
					//diff = undriven - ref;
					// Any of this mean it's not a valid sample to consider for zero crossing
					if(high > low || abs(diff) > 150)
						continue;
					running_avg = 0.7 * running_avg + 0.3 * diff; //(undriven - ref);
					//diff = running_avg;
				} else {
					diff = MID_POINT - undriven - state_offset[curr_state];
					
					big_val = (ref - undriven);
					// If either of these true not a good sample to consider
					if(high < low || diff > 150) 
						continue;					
					
				}
				
				if(diff < 0) {
					pre_count++;
					// Keep setting so we store the time of zero crossing
					below_time = enter_time - dT * (DOWNSAMPLING - i);
				}
				if(diff >0 && pre_count > MIN_PRE_COUNT) {
					post_count++;
				}
				if(diff > 0 && post_count >= MIN_POST_COUNT) {
					//below_time = enter_time;// - dT * (DOWNSAMPLING - i);
					commutation_detected = true;
					break;
				}				
			}
			break;
		default:
			//PIOS_ESC_Off();
			break;
	}
	
	
	
	if(commutation_detected) {		
		
		if(message.read == false)
			PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG,(uint8_t *) "Oh fuck\n", 8);
		
		message.state = curr_state;
		message.time = below_time;		
		message.read = false;
		process_message(&message);
	}
	
//#define DUMP_ADC
#ifdef DUMP_ADC
	static uint8_t count = 0;
	uint8_t buf[4 + DOWNSAMPLING * 4 * 2];
	buf[0] = 0x00; // syncing bytes
	buf[1] = 0xff;
	buf[2] = 0xc5;
	buf[3] = count++;
	memcpy(&buf[4], (uint8_t *) raw_buf, DOWNSAMPLING * 4 * 2);
	
	PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, buf, sizeof(buf));
#endif
	
	
		
//#define DUMP_UNDRIVEN
#ifdef DUMP_UNDRIVEN
	for(int i = 0; i < DOWNSAMPLING; i++) {
		uint32_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
		uint32_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
		
		if((high > low) && (high < (low + DRIVING_THRESHOLD))) {
			uint32_t sample = PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin;
			PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, (uint8_t *) &raw_buf[sample], 2);
		}
	}
#endif
	
//#define DUMP_DIFF_DRIVEN
#ifdef DUMP_DIFF_DRIVEN
	// This seems to be some decent processing to get downward slope commutation detection
	for(int i = 0; i < DOWNSAMPLING; i++) {
		int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
		int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
		int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
		bool bad = false;
		int16_t ref = (high + low) / 2; 
		if(high < low || abs(undriven - ref) > 100) {
			bad = true;
		}
		if(!bad && curr_state == ESC_STATE_AB) {
			uint16_t buf[2];
			buf[0] = undriven;
			buf[1] = ref;
			PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, (uint8_t*) &buf, sizeof(buf));
		}

	}
#endif
	
//#define DUMP_DRIVEN
#ifdef DUMP_DRIVEN
	static uint8_t count = 0;
	uint8_t buf[4 + DOWNSAMPLING * 2];
	buf[0] = 0x00; // syncing bytes
	buf[1] = 0xff;
	buf[2] = 0xc3;
	buf[3] = count++;
	for(int i = 0; i < DOWNSAMPLING; i++) {
		uint32_t sample = PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin;
		buf[4 + i * 2] = raw_buf[sample] & 0xff;
		buf[4 + i * 2 + 1] = (raw_buf[sample] >> 8) & 0xff;
	}
	PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, buf, sizeof(buf));			
#endif

adc_done:
	adc_irq_status = false;

}
/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still 
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


