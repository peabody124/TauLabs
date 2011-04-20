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

//TODO: Slave two timers together so in phase

/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

void adc_callback(float * buffer);

#define DOWNSAMPLING 6

// A message from the ADC to say a zero crossing was detected
struct zerocrossing_message {
	enum pios_esc_state state;
	uint16_t time;
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

// TODO: A tim4 interrupt that actually implements the commutation on a regular schedule

#define COMMUTATIONS_PER_ROT (7*6)
float initial_startup_speed = 150;
float final_startup_speed = 1500;
float current_speed;
bool closed_loop = false;
int32_t desired_closed_delay = 800;
int32_t desired_delay = 500;
volatile int32_t delay = 5000;
volatile uint16_t next;
volatile uint16_t swap_time;
volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin; 
volatile bool pos;
volatile float dc = 0.13;

volatile bool commutation_detected = false;
volatile uint16_t checks = 0;
uint16_t consecutive_nondetects = 0;

uint8_t dT = 6; // 6 uS per sample at 160k
float rate = 0;
uint16_t skipped_delay = 1050;

void process_message(struct zerocrossing_message * msg);

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

	// This pull up all the ADC voltages so the BEMF when at -0.7V
	// is still positive
	PIOS_GPIO_Enable(0);
	PIOS_GPIO_Off(0);
		
	PIOS_LED_Off(LED1);
	PIOS_LED_On(LED2);
	PIOS_LED_On(LED3);
	
	int32_t count = 0;
	
	next = PIOS_DELAY_GetuS() + delay;

	PIOS_ESC_SetDutyCycle(dc);
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);
	PIOS_ESC_Arm();
	
	
	while(1) {
		//process_message(&message);
			
		if(PIOS_DELAY_DiffuS(next) > 0) {
			//PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", (next - swap_time), (uint32_t) (dc * 10000));
			swap_time = PIOS_DELAY_GetuS();
			PIOS_ESC_NextState();

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
					if(fail_count++ > 50) {
						PIOS_ESC_Off();
						break;
					}
				}
//				else
//					PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", PIOS_ESC_GetState(), 100 * zerocrossing_stats.latency[PIOS_ESC_GetState()] / zerocrossing_stats.smoothed_interval);
				closed_loop_updated = false;

				if(consecutive_nondetects > 10)
					PIOS_ESC_Off();

				// Update duty cycle and such 
				next += skipped_delay;
//				PIOS_ESC_SetDutyCycle(0.18);
//				dc += 0.000001 * (rate - desired_closed_delay);
//				if(dc > 0.08 && dc < 0.20)
//				PIOS_ESC_SetDutyCycle(dc); 
/*				dc = 0.12 + (dc - 0.12) * 0.999;
				PIOS_ESC_SetDutyCycle(dc); */
			} else {
				// Turn err light off
				PIOS_LED_Off(LED2);

				//TODO: Run startup state machine

				if(current_speed < final_startup_speed) 
					current_speed+=0.4;
				if(current_speed >= final_startup_speed) {
					dc = 0.1 + (dc - 0.1) * 0.9999;
					PIOS_ESC_SetDutyCycle(dc);
				}
				
				delay = 1e6 * 60 / (current_speed * COMMUTATIONS_PER_ROT);
				next += delay;
				
			}			

			count++;									
			
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
		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u %u\n", msg->state, skipped_delay, zerocrossing_stats.smoothed_interval);

	// If decent interval use it to update estimate of speed
	if(!skipped && !prev_skipped && (zerocrossing_stats.interval < 2000))
		zerocrossing_stats.smoothed_interval = 0.9 * zerocrossing_stats.smoothed_interval + 0.1 * zerocrossing_stats.interval;

//		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", zerocrossing_stats.smoothed_interval, zerocrossing_stats.interval);
	if(zerocrossing_stats.consecutive_detected > 100) 
		closed_loop = true;
							
	if(closed_loop) {
		// TOOD: This logic shouldn't stay here
		closed_loop_updated = true;
		skipped_delay = delay;
			//uint16_t pre_next = next;
		next = msg->time + zerocrossing_stats.smoothed_interval * 0.40;
		skipped_delay = zerocrossing_stats.smoothed_interval * 20;
			//PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u %u %u %u\n", swap_time, pre_next, next, zerocrossing_stats.smoothed_interval, zerocrossing_stats.interval);
	}

	
	prev_skipped = skipped;
	prev_state = msg->state;
}

// When driving both legs to ground mid point is 580
#define MID_POINT 580
#define DIODE_LOW 460
#define MIN_PRE_COUNT  2
#define MIN_POST_COUNT 2
#define DEMAG_BLANKING 20
#define UPSLOPE_OVERSHOOT 80
#define DOWNSLOPE_OVERSHOOT 60

void adc_callback(float * buffer) 
{
	int16_t * raw_buf = PIOS_ADC_GetRawBuffer();
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	enum pios_esc_state curr_state;	
	static uint16_t below_time;
	static uint16_t pre_count = 0;
	static uint16_t post_count = 0;
	static int16_t running_avg;

	// Commutation detection, assuming mode is ESC_MODE_LOW_ON_PWM_BOTH
	uint16_t enter_time = PIOS_DELAY_GetuS();
	// Process ADC from undriven leg
	// TODO: Make these variables be updated when a state transition occurs

	// Wait for blanking
	if(((uint16_t) (enter_time - swap_time)) < DEMAG_BLANKING)
		return;
	
	if(commutation_detected)
		return;
	
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
			for(int i = 0; i < DOWNSAMPLING; i++) {
				int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
				int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
				int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
				//					int16_t ref = (high + MID_POINT) / 2; 
				int16_t diff;
				
				// For now only processing the low phase of the duty cycle
				if(high > 3000)
					continue; 
				
				if(pos) {						
					diff = undriven - MID_POINT - DOWNSLOPE_OVERSHOOT; 
					
					// Any of this mean it's not a valid sample to consider for zero crossing
					if(high > low || high > 1000 || abs(diff) > 120)
						continue;
					running_avg = 0.7 * running_avg + 0.3 * diff; //(undriven - ref);
					//						diff = running_avg;
				} else {
					diff = MID_POINT - undriven - UPSLOPE_OVERSHOOT;
					
					// If either of these true not a good sample to consider
					if(high < low || abs(diff) > 120) 
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
	
}
/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still 
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


