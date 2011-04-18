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


/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

void adc_callback(float * buffer);

#define DOWNSAMPLING 12

// TODO: A tim4 interrupt that actually implements the commutation on a regular schedule

int32_t desired_delay = 1500;
volatile int32_t delay = 4000;
volatile uint16_t next;
volatile uint16_t swap_time;
volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin; 
volatile bool pos;
volatile float dc = 0.1;
uint8_t dT = 10; // 10 uS per sample at 100k

/**
 * @brief ESC Main function
 */
int main()
{
	PIOS_Board_Init();
	
	PIOS_ADC_Config(DOWNSAMPLING);
	PIOS_ADC_SetCallback(adc_callback); 

	// This pull up all the ADC voltages so the BEMF when at -0.7V
	// is still positive
	PIOS_GPIO_Enable(0);
	PIOS_GPIO_Off(0);
	
	uint8_t closed_loop = 0;
	uint8_t commutation_detected = 0;
	
	PIOS_LED_Off(LED1);
	PIOS_LED_On(LED2);
	PIOS_LED_On(LED3);
	
	int32_t count = 0;
	
	next = PIOS_DELAY_GetuS() + delay;

	PIOS_ESC_SetDutyCycle(dc);
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);
	PIOS_ESC_Arm();
	
	
	while(1) {
		//Process analog data, detect commutation
		commutation_detected = 0;
		
		if(closed_loop) {
			// Update duty cycle and such 
		} else {
			// Run startup state machine
		}
		
			
		if(PIOS_DELAY_DiffuS(next) > 0) {
			
			next += delay;
		
			PIOS_LED_Toggle(LED3);

			count++;
			
			if(delay > desired_delay)
				delay--;
						
			PIOS_ESC_NextState();
			swap_time = PIOS_DELAY_GetuS();
			
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
		
	}
	return 0;
}	

// When driving both legs to ground mid point is 580
#define MID_POINT 580
#define DIODE_LOW 460
#define DRIVING_THRESHOLD 400
#define COMMUTATION_THRESHOLD 100
#define RUNNING_FILTER_SIZE 3
#define MIN_PRE_COUNT  3
#define MIN_POST_COUNT 3

void adc_callback(float * buffer) 
{
	int16_t * raw_buf = PIOS_ADC_GetRawBuffer();
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	enum pios_esc_state curr_state;	
	static bool swapped = false;
	static bool detected = false;
	static uint16_t last_time;
	static uint16_t this_time;
	static uint16_t below_time;
	static bool outputted = false;
	static bool active = false;
	static uint16_t detection_count = 0;
	static uint16_t pre_count = 0;
	static uint16_t post_count = 0;
	
	static bool started = false;

	// Commutation detection, assuming mode is ESC_MODE_LOW_ON_PWM_BOTH
	curr_state = PIOS_ESC_GetState();
	uint16_t enter_time = PIOS_DELAY_GetuS();
	// Process ADC from undriven leg
	// TODO: Make these variables be updated when a state transition occurs

	
	if(curr_state != prev_state) {
		prev_state = curr_state;
		swapped = true;
		detected = false;
		outputted = false;
		started = false;
		swap_time = PIOS_DELAY_GetuS();

		pre_count = 0;
		post_count = 0;
	}

	if(swapped && !detected) { // Once detected don't need to process till next state
		switch(PIOS_ESC_GetMode()) {
			case ESC_MODE_LOW_ON_PWM_HIGH:
				// Doesn't work quite right yet
				for(int i = 0; i < DOWNSAMPLING; i++) {
					int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
					int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
					int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
					int16_t diff = undriven - low; 
					
					if(pos) {
						// Any of this mean it's not a valid sample to consider for zero crossing
						if((high > low) ||
						   (undriven > 1000) || 
						   (high > 1000)  ||
						   (low > 1000)  ||
						   (abs(diff) > 60))
							continue;
																		
						if(diff < 0) {
							pre_count++;
							// Keep setting so we store the time of zero crossing
							below_time = enter_time - dT * (DOWNSAMPLING - i);
						}
						if(diff >0 && pre_count > MIN_PRE_COUNT) {
							post_count++;
						}
						if(diff > 0 && post_count >= MIN_POST_COUNT) {
							if(!detected)
								detection_count++;
							detected = true;
						}
					}
				}
				break;
			case ESC_MODE_LOW_ON_PWM_BOTH:
				// Doesn't work quite right yet
				for(int i = 0; i < DOWNSAMPLING; i++) {
					uint32_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
					uint32_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
					if(abs(high - MID_POINT) < DRIVING_THRESHOLD && abs(low - MID_POINT) < DRIVING_THRESHOLD) {
						// If this sample is from when both driving sides pulled low
						uint32_t sample = PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin;
						if(abs(raw_buf[sample] - MID_POINT) > COMMUTATION_THRESHOLD) {
							if(!detected)
								detection_count++;
							detected = true;
						}
					}
				}
				break;
			default:
				PIOS_ESC_Off();
		}
	}
	
	
	if(detection_count > 100) 
		active = true;

	if(detected && !outputted && active) {
		last_time = this_time;
		this_time = below_time;
		uint16_t wait = this_time - swap_time;
		
		static float running_avg = 1000;
		running_avg = running_avg * 0.65 + (float) wait * 0.35;
				
		dc += (0.000001 * (2 * running_avg - desired_delay));
		if(dc < 0.20 && dc  > 0.05)
			PIOS_ESC_SetDutyCycle(dc);
		else 
			dc = 0.05;
/*
		if(detection_count > 1000) {
			if(delay > 2 * running_avg)
				delay --;
			if(delay < 2 * running_avg)
				delay++;
		} */
				
		
		// If tracking closely refine next commutation time
//		if(abs(delay - 2 * running_avg) < 100)
//			next = this_time + delay / 2;
			
		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u %u %u\n", delay, (uint16_t) this_time - swap_time, (uint16_t) next - this_time, (int32_t) running_avg);
		
//		next = PIOS_DELAY_GetuS() + (this_time - state_swap_time) * 1.0;
//		next = PIOS_DELAY_GetuS() + delay;
		

		outputted = true;
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
	// This seems to be some decent processing to get upward slope commutation detection
	static uint8_t count = 0;
	uint8_t buf[4 + DOWNSAMPLING * 2];
	buf[0] = 0x00; // syncing bytes
	buf[1] = 0xff;
	buf[2] = 0xc7;
	buf[3] = count++;
	for(int i = 0; i < DOWNSAMPLING; i++) {
		//int16_t high = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + high_pin];
		int16_t low = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + low_pin];
		int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin];
//		int16_t ref = (high + low) / 2; 
		int16_t diff = undriven - low; 
		buf[4 + i * 2] = diff & 0xff;
		buf[4 + i * 2 + 1] = (diff >> 8) & 0xff;
/*		if(high < low) {
			buf[4 + i * 2] = 0;
			buf[4 + i * 2 + 1] = 0;
		} 
		if(undriven > 1000) {
			buf[4 + i * 2] = 1;
			buf[4 + i * 2 + 1] = 0;
		}
		if(high > 1000) {
			buf[4 + i * 2] = 2;
			buf[4 + i * 2 + 1] = 0;
		} 
		if(low > 1000) {
			buf[4 + i * 2] = 3;
			buf[4 + i * 2 + 1] = 0;
		}
		if(abs(diff) > 50) {
			buf[4 + i * 2] = 4;
			buf[4 + i * 2 + 1] = 0;
		} */
	}
	if(curr_state == ESC_STATE_AB)
	PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, buf, sizeof(buf));			
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


