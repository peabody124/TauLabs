/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup PIOS_DELAY Delay Functions
 * @brief PiOS Delay functionality
 * @{
 *
 * @file       pios_esc.h
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

#ifndef PIOS_ESC_H
#define PIOS_ESC_H

typedef enum esc_bridge_mode {LOW_ON_PWM_HIGH = 0, LOW_ON_PWM_LOW = 1, LOW_ON_PWM_BOTH = 2, HIGH_ON_PWM_LOW = 3, HIGH_ON_PWM_HIGH = 4, HIGH_ON_PWM_BOTH = 5};
typedef enum esc_state {ESC_STATE_AB = 0, ESC_STATE_AC = 1, ESC_STATE_BA = 2, ESC_STATE_BC = 3, ESC_STATE_CA = 4, ESC_STATE_CB = 5};
void PIOS_ESC_SetState(esc_state new_state, uint16_t duty_cycle, esc_bridge_mode mode = 0);

#endif /* PIOS_DELAY_H */

/**
 * @}
 * @}
 */

