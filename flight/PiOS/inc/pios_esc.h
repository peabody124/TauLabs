/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup PIOS_ESC ESC Functions
 * @brief PiOS ESC functionality
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

//TODO: Add ID to support multiple ESC outputs
void PIOS_ESC_Off();
void PIOS_ESC_Arm();
void PIOS_ESC_NextState();
void PIOS_ESC_SetDutyCycle(uint16_t duty_cycle);
void PIOS_ESC_SetState(uint8_t new_state);

#endif /* PIOS_DELAY_H */

/**
 * @}
 * @}
 */

