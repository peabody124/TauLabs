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
};

// Circular buffer can hold two commutation messages
uint8_t message_buffer_space[sizeof(struct zerocrossing_message) * 2];
t_fifo_buffer message_buffer;

// TODO: A tim4 interrupt that actually implements the commutation on a regular schedule

int32_t desired_delay = 950;
volatile int32_t delay = 3000;
volatile uint16_t next;
volatile uint16_t swap_time;
volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin; 
volatile bool pos;
volatile float dc = 0.11;
float desired_dc = 0.125;
uint8_t dT = 10; // 10 uS per sample at 100k
float rate = 0;
uint16_t skipped_delay = 950;
/**
 * @brief ESC Main function
 */
int main()
{
	// Must init the buffer before the ADC
	fifoBuf_init(&message_buffer, message_buffer_space, sizeof(message_buffer_space));
	
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
				
		if(fifoBuf_getUsed(&message_buffer) >= sizeof(struct zerocrossing_message)) {
			// There is a message waiting
//			static int dropped = 0;
			static int detected = 0;
			bool skipped = false;
			static bool prev_skipped = false;
			static uint32_t skipped_count = 0;
			static uint32_t consecutive_skipped = 0;
			static uint16_t last_time;
			static enum pios_esc_state prev_state;
			struct zerocrossing_message message;
			switch(message.state) {
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
			prev_state = message.state;
			
			if(!skipped)
				consecutive_skipped = 0;
			else 
				consecutive_skipped++;
			
			if(consecutive_skipped > 5) {
				detected = 0;
				closed_loop = false;
			}
				
			fifoBuf_getData(&message_buffer, &message, sizeof(message));
			uint16_t interval = message.time - last_time;
			PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u %u %u\n", message.state, interval, skipped_count, (uint32_t) rate);
			
			last_time = message.time;
			
			uint16_t buf[4];
			buf[0] = 0x03;
			buf[1] = message.state;
			buf[2] = interval;
			buf[3] = message.time - swap_time;
//			PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, (uint8_t *) &buf[0], sizeof(buf));
			
			if( ! skipped && ! prev_skipped && (interval < 1000)) {				
				rate = 0.8 * rate + 0.2* interval;
				if(detected > 100) {					
					closed_loop = true;
					
					uint16_t new_next = message.time + 210; //100 + rate * 0.2;
//					if(abs((int16_t) (new_next - next)) < 200) {
						next = new_next;
//					}
				} 
				if(detected > 3000)
					skipped_delay = 950; //rate + 100;

				detected++;
				
				if(rate < 700)
					dc += 0.0001;
			/*
			dc += (0.00000001 * (rate - 2 * delay));
			if(dc < 0.20 && dc  > 0.05)
				PIOS_ESC_SetDutyCycle(dc);
			else 
				dc = 0.05;    */
			} else { 
				skipped_count++;
			} 	
			
			prev_skipped = skipped;
			
		}
		
			
		if(PIOS_DELAY_DiffuS(next) > 0) {
			PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u\n", (next - swap_time));
			swap_time = PIOS_DELAY_GetuS();
			PIOS_ESC_NextState();

			if(closed_loop) {
				// Update duty cycle and such 
				next += skipped_delay;
				
				if(desired_dc > dc) {
					float diff = desired_dc - dc;
					diff *= 0.001;
					if(diff > 0.000001)
						diff = 0.000001;
//					dc += diff;
//					PIOS_ESC_SetDutyCycle(dc);
				}
			} else {
				//TODO: Run startup state machine
				next += delay;
				
				if(delay > desired_delay)
					delay -= 2;				
			}			
		
			PIOS_LED_Toggle(LED3);

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
		}			
		
	}
	return 0;
}	

// When driving both legs to ground mid point is 580
#define MID_POINT 580
#define DIODE_LOW 460
#define MIN_PRE_COUNT  2
#define MIN_POST_COUNT 2
#define DEMAG_BLANKING 50

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

	// Wait for blanking
	if(((uint16_t) (enter_time - swap_time)) < DEMAG_BLANKING)
		return;
	
	// If detected this commutation don't bother here
	// TODO:  disable IRQ for efficiency
	if(curr_state == prev_state && detected)
		return;

	
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
//					int16_t ref = (high + MID_POINT) / 2; 
					int16_t diff;
					
					// For now only processing the low phase of the duty cycle
					if(high > 3000)
						continue; 
					
					if(pos) {						
						diff = undriven - MID_POINT; 

						// Any of this mean it's not a valid sample to consider for zero crossing
						if(high > low || high > 1000 || abs(diff) > 120)
							continue;
					} else {
						diff = MID_POINT - undriven;
						
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
						if(!detected)
							detection_count++;
						detected = true;
					}
					
				}
				break;
			default:
				//PIOS_ESC_Off();
				break;
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
				
		struct zerocrossing_message message = {
			.state = curr_state,
			.time = below_time
		};
		fifoBuf_putData(&message_buffer, &message, sizeof(message));
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
			
//		PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u %u %u\n", delay, (uint16_t) this_time - swap_time, (uint16_t) next - this_time, (int32_t) running_avg);
		
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


