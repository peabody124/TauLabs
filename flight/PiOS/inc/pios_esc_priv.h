/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup PIOS_ESC Priviledged definitions
 * @brief PiOS ESC functionality
 * @{
 *
 * @file       pios_esc_priviledged.h
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * 	       James Cotton
 * @brief      ESC functions header 
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

#ifndef PIOS_ESC_PRIV_H
#define PIOS_ESC_PRIV_H

#define ESC_MODE_LOW_ON_PWM_HIGH  0
#define ESC_MODE_LOW_ON_PWM_LOW   1
#define ESC_MODE_LOW_ON_PWM_BOTH  2
#define ESC_MODE_HIGH_ON_PWM_LOW  3
#define ESC_MODE_HIGH_ON_PWM_HIGH 4
#define ESC_MODE_HIGH_ON_PWM_BOTH 5

#define ESC_STATE_AB              0
#define ESC_STATE_AC              1
#define ESC_STATE_BA              2
#define ESC_STATE_BC              3
#define ESC_STATE_CA              4
#define ESC_STATE_CB              5

//! Macros to determine if 100% driven leg is high or low
#define ESC_LOW_ON(mode)  (mode == ESC_MODE_LOW_ON_PWM_HIGH || mode == ESC_MODE_LOW_ON_PWM_LOW || mode == ESC_MODE_LOW_ON_PWM_BOTH)
#define ESC_HIGH_ON(mode)  (mode == ESC_MODE_HIGH_ON_PWM_LOW || mode == ESC_MODE_HIGH_ON_PWM_HIGH || mode == ESC_MODE_HIGH_ON_PWM_BOTH)

//! Macros to determine if the high or low side (or both) should be PWM'd
#define ESC_PWM_HIGH(mode) (mode == ESC_MODE_LOW_ON_PWM_HIGH || mode == ESC_MODE_HIGH_ON_PWM_HIGH \
|| mode == ESC_MODE_LOW_ON_PWM_BOTH || mode == ESC_MODE_HIGH_ON_PWM_BOTH)
#define ESC_PWM_LOW(mode) (mode == ESC_MODE_LOW_ON_PWM_LOW || mode == ESC_MODE_HIGH_ON_PWM_LOW \
|| mode == ESC_MODE_LOW_ON_PWM_BOTH || mode == ESC_MODE_HIGH_ON_PWM_BOTH)

//! Macros to determine if a leg is just driven on (100% duty cycle)
#define ESC_A_STEADY(state) (state == ESC_STATE_AB || state == ESC_STATE_AC)
#define ESC_B_STEADY(state) (state == ESC_STATE_BA || state == ESC_STATE_BC)
#define ESC_C_STEADY(state) (state == ESC_STATE_CA || state == ESC_STATE_CB)

//! Macros to determine if a leg is driven pwm
#define ESC_A_PWM(state) (state == ESC_STATE_BA || state == ESC_STATE_CA)
#define ESC_B_PWM(state) (state == ESC_STATE_AB || state == ESC_STATE_CB)
#define ESC_C_PWM(state) (state == ESC_STATE_AC || state == ESC_STATE_BC)

//! Macros to determine if a leg is not driven
#define ESC_A_OFF(state) (state == ESC_STATE_BC || state == ESC_STATE_CB)
#define ESC_B_OFF(state) (state == ESC_STATE_AC || state == ESC_STATE_CA)
#define ESC_C_OFF(state) (state == ESC_STATE_AB || state == ESC_STATE_BA)

struct pios_esc_cfg {
	// PWM timer configuration
	// Pin mappings
};

struct pios_esc_dev {
	const struct pios_esc_cfg *const cfg;	
	uint8_t state;
	uint8_t mode;
	uint8_t armed;
	uint16_t duty_cycle;
};

extern struct pios_esc_dev pios_esc_devs[];
extern uint8_t pios_esc_num_devices;

#endif