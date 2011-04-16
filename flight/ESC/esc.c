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

#define DOWNSAMPLING 5

/**
 * @brief ESC Main function
 */
int main()
{
	PIOS_Board_Init();
	
	PIOS_ADC_Config(DOWNSAMPLING);
	PIOS_ADC_SetCallback(adc_callback); 

	PIOS_ESC_SetDutyCycle(0.12);
	//PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_BOTH);
	PIOS_ESC_Arm();

	uint8_t closed_loop = 0;
	uint8_t step = 1;  // for testing leave at 1
	uint8_t commutation_detected = 0;
	
	PIOS_LED_Off(LED1);
	PIOS_LED_On(LED2);
	PIOS_LED_On(LED3);
	
	int32_t delay = 5000;
	uint16_t next = PIOS_DELAY_GetuS() + delay;
	
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

			if( delay > 980)  {
				delay --;
			} else {
				PIOS_ESC_Off();
				while(1);
			}
			
			if(step) 
				PIOS_ESC_NextState();
		}			
		
	}
	return 0;
}	

void adc_callback(float * buffer) 
{
	static uint8_t count = 0;
	uint8_t buf[4 + DOWNSAMPLING * 4 * 2];
	buf[0] = 0x00; // syncing bytes
	buf[1] = 0xff;
	buf[2] = 0xc3;
	buf[3] = count++;
	memcpy((uint8_t *) PIOS_ADC_GetRawBuffer(), &buf[4], DOWNSAMPLING * 4 * 2);
	
	PIOS_COM_SendBufferNonBlocking(PIOS_COM_DEBUG, buf, sizeof(buf));
}
/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still 
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


