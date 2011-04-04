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

#if defined(PIOS_INCLUDE_ESC)

/*
 
 struct pios_esc_cfg {
 // PWM timer configuration
 // Pin mappings
 };
 
 struct pios_esc_dev {
 const struct pios_esc_cfg *const cfg;	
 uint8_t state;
 uint8_t mode;
 uint8_t armed;
 };
 
 extern struct pios_esc_dev pios_esc_devs[];
 extern uint8_t pios_esc_num_devices;
*/


/**
 * @brief Initialize the driver outputs for PWM mode.  Configure base timer rate.
 */
void PIOS_Init()
{
}

/**
 * @brief Function for immediately shutting off the bridge (not a braking mode) in case of faults
 * also disarms
 */
void PIOS_ESC_Off()
{
	//TOOD: Disarm
}

/**
 * @brief Arm the ESC.  Required for any driver outputs to go high.
 */
void PIOS_ESC_Arm()
{
}

/**
 * @brief Take a step through the state chart of modes
 */
void PIOS_ESC_NextState() 
{
}

/**
 * @brief Sets the duty cycle of all PWM outputs
 */
void PIOS_ESC_SetDutyCycle(uint16_t duty_cycle)
{
}

/**
 * @brief Set the ESC to the new state
 * @param[in] new_state Determines which pair of arms are energized
 * @param[in] duty_cycle The amount of time to energize for
 * @param[in] mode Allows selecting whether the high or low side is switched 
 * as well as active braking 
 */ 
void PIOS_ESC_SetState(esc_state new_state, uint16_t duty_cycle, esc_bridge_mode mode = 0) 
{
}

#endif