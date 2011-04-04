/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup PIOS_ESC ESC Functions
 * @brief PiOS ESC functionality
 * @{
 *
 * @file       pios_esc.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * 	       James Cotton
 * @brief      ESC Functions 
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

/* Project Includes */
#include "pios.h"
#include "pios_esc_priv.h"
#include <stdbool.h>

#if defined(PIOS_INCLUDE_ESC)

//! Private functions
static void PIOS_ESC_SetMode(uint8_t mode);


/**
 * @brief Initialize the driver outputs for PWM mode.  Configure base timer rate.
 */
void PIOS_Init()
{
	PIOS_ESC_Off();
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);
}

/**
 * @brief Function for immediately shutting off the bridge (not a braking mode) in case of faults
 * also disarms
 */
void PIOS_ESC_Off()
{
	pios_esc_devs[0].armed = false;
	PIOS_ESC_SetDutyCycle(0);
}

/**
 * @brief Arm the ESC.  Required for any driver outputs to go high.
 */
void PIOS_ESC_Arm()
{
	pios_esc_devs[0].armed = true;
}

/**
 * @brief Take a step through the state chart of modes
 */
void PIOS_ESC_NextState() 
{
	uint8_t new_state;
	
	// Simple FSM
	switch(pios_esc_devs[0].state) {
		case ESC_STATE_AB:
			new_state = ESC_STATE_CB;
			break;
		case ESC_STATE_AC:
			new_state = ESC_STATE_AB;
			break;
		case ESC_STATE_BA:
			new_state = ESC_STATE_BC;
			break;
		case ESC_STATE_BC:
			new_state = ESC_STATE_AC;
			break;
		case ESC_STATE_CA:
			new_state = ESC_STATE_BA;
			break;
		case ESC_STATE_CB:
			new_state = ESC_STATE_CA;
			break;
		default: // should not happen
			new_state = ESC_STATE_AB;
	}
			
	PIOS_ESC_SetState(new_state);
}

/**
 * @brief Sets the duty cycle of all PWM outputs
 */
void PIOS_ESC_SetDutyCycle(uint16_t duty_cycle)
{
	pios_esc_devs[0].duty_cycle = duty_cycle;
	// TODO: Update the outputs accordingly
}

/**
 * @brief Sets the duty cycle of all PWM outputs
 */
static void PIOS_ESC_SetMode(uint8_t mode)
{
	if(mode <= ESC_MODE_HIGH_ON_PWM_BOTH)
		pios_esc_devs[0].mode = mode;
}

/**
 * @brief Set the ESC to the new state
 * @param[in] new_state Determines which pair of arms are energized
 * @param[in] duty_cycle The amount of time to energize for
 * @param[in] mode Allows selecting whether the high or low side is switched 
 * as well as active braking 
 */ 
void PIOS_ESC_SetState(uint8_t new_state)
{
}

#endif