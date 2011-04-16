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
#include "pios_esc.h"

enum pios_esc_state {ESC_STATE_AB = 0,
	ESC_STATE_AC,
	ESC_STATE_BA,
	ESC_STATE_BC,
	ESC_STATE_CA,
	ESC_STATE_CB};

enum pios_esc_gate_mode {ESC_GATE_MODE_OFF = 0,
	ESC_GATE_MODE_ON,
	ESC_GATE_MODE_PWM,
	ESC_GATE_MODE_PWM_INVERT};

#define ESC_PHASE_A               1
#define ESC_PHASE_B               2
#define ESC_PHASE_C               3

struct pios_gate_channel {
	TIM_TypeDef * timer;
	GPIO_TypeDef * port;
	uint8_t channel;
	uint16_t pin;
};

struct pios_esc_cfg {
	TIM_TimeBaseInitTypeDef tim_base_init;
	TIM_OCInitTypeDef tim_oc_init;
	GPIO_InitTypeDef gpio_init;
	uint32_t remap;
	uint16_t pwm_base_rate;
	struct pios_gate_channel phase_a_minus;
	struct pios_gate_channel phase_a_plus;
	struct pios_gate_channel phase_b_minus;
	struct pios_gate_channel phase_b_plus;
	struct pios_gate_channel phase_c_minus;
	struct pios_gate_channel phase_c_plus;
};

struct pios_esc_dev {
	const struct pios_esc_cfg * cfg;	
	enum pios_esc_state state;
	enum pios_esc_mode mode;
	volatile uint8_t armed;
	volatile uint16_t duty_cycle;
	volatile enum pios_esc_gate_mode phase_a_minus;
	volatile enum pios_esc_gate_mode phase_a_plus;
	volatile enum pios_esc_gate_mode phase_b_minus;
	volatile enum pios_esc_gate_mode phase_b_plus;
	volatile enum pios_esc_gate_mode phase_c_minus;
	volatile enum pios_esc_gate_mode phase_c_plus;
};

extern struct pios_esc_dev pios_esc_dev;

void PIOS_ESC_Init(const struct pios_esc_cfg * cfg);

#endif