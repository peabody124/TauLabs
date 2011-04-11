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

/**
 * @brief ESC Main function
 */
int main()
{
	PIOS_Board_Init();
	
	PIOS_ESC_SetDutyCycle(10);
	PIOS_ESC_Arm();
	
	uint8_t closed_loop = 0;
	uint8_t step = 1;  // for testing leave at 1
	uint8_t commutation_detected = 0;
	
	PIOS_LED_On(1);
	while(1) {
		
		//Process analog data, detect commutation
		commutation_detected = 0;
		
		if(closed_loop) {
			// Update duty cycle and such 
		} else {
			// Run startup state machine
		}
		
			
		if(step) 
			PIOS_ESC_NextState();

		PIOS_DELAY_WaituS(1000);
	}
	return 0;
}	

/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still 
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


