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
static void PIOS_ESC_UpdateOutputs();
static void PIOS_ESC_UpdateChannel(const struct pios_gate_channel * gate, uint8_t gate_mode);

struct pios_esc_dev pios_esc_dev;

#define LAST_STATE ESC_STATE_CB
#define LAST_MODE  ESC_MODE_HIGH_ON_PWM_BOTH

// TOOD: Add supervisor that times out and disarms

/**
 * @brief Initialize the driver outputs for PWM mode.  Configure base timer rate.
 */
void PIOS_ESC_Init(const struct pios_esc_cfg * cfg)
{
	pios_esc_dev.cfg = cfg;
	PIOS_ESC_Off();
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);
	
	
	for (uint8_t i = 0; i < 6; i++) {
		GPIO_InitTypeDef GPIO_InitStructure = cfg->gpio_init;
		TIM_TimeBaseInitTypeDef TIM_TimeBaseStructure = cfg->tim_base_init;
		TIM_OCInitTypeDef TIM_OCInitStructure = cfg->tim_oc_init;
				
		struct pios_gate_channel channel;
		switch(i) {
			case 0:
				channel = cfg->phase_a_minus;
				break;
			case 1:
				channel = cfg->phase_a_plus;
				break;
			case 2:
				channel = cfg->phase_b_minus;
				break;
			case 3:
				channel = cfg->phase_b_plus;
				break;
			case 4:
				channel = cfg->phase_c_minus;
				break;
			case 5:
				channel = cfg->phase_c_plus;
				break;
			default:
				// Serious error
				PIOS_ESC_Off();
				return;
		}
		
		/* Enable appropriate clock to timer module */
		switch((int32_t) channel.timer) {
			case (int32_t)TIM1:
				RCC_APB2PeriphClockCmd(RCC_APB2Periph_TIM1, ENABLE);
				break;
			case (int32_t)TIM2:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM2, ENABLE);
				break;
			case (int32_t)TIM3:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM3, ENABLE);
				break;
			case (int32_t)TIM4:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM4, ENABLE);
				break;
			case (int32_t)TIM5:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM5, ENABLE);
				break;
			case (int32_t)TIM6:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM6, ENABLE);
				break;
			case (int32_t)TIM7:
				RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM7, ENABLE);
				break;
			case (int32_t)TIM8:
				RCC_APB2PeriphClockCmd(RCC_APB2Periph_TIM8, ENABLE);
				break;
		}
		
		/* Enable GPIO */
		GPIO_InitStructure.GPIO_Pin = channel.pin;
		GPIO_Init(channel.port, &GPIO_InitStructure);
		
		/* Enable time base */
		TIM_TimeBaseInit(channel.timer,  &TIM_TimeBaseStructure);
		
		channel.timer->PSC = cfg->tim_base_init.TIM_Prescaler;
		
		/* Set up for output compare function */
		switch(channel.channel) {
			case TIM_Channel_1:
				TIM_OC1Init(channel.timer, &TIM_OCInitStructure);
				TIM_OC1PreloadConfig(channel.timer, TIM_OCPreload_Enable);
				break;
			case TIM_Channel_2:
				TIM_OC2Init(channel.timer, &TIM_OCInitStructure);
				TIM_OC2PreloadConfig(channel.timer, TIM_OCPreload_Enable);
				break;
			case TIM_Channel_3:
				TIM_OC3Init(channel.timer, &TIM_OCInitStructure);
				TIM_OC3PreloadConfig(channel.timer, TIM_OCPreload_Enable);
				break;
			case TIM_Channel_4:
				TIM_OC4Init(channel.timer, &TIM_OCInitStructure);
				TIM_OC4PreloadConfig(channel.timer, TIM_OCPreload_Enable);
				break;
		}
		
		RCC_APB2PeriphClockCmd(RCC_APB2Periph_GPIOA, ENABLE);
		TIM_ARRPreloadConfig(channel.timer, ENABLE);
		TIM_CtrlPWMOutputs(channel.timer, ENABLE);
		TIM_Cmd(channel.timer, ENABLE);		
		
	}	
	
	if(cfg->remap) {
		/* Warning, I don't think this will work for multiple remaps at once */
		GPIO_PinRemapConfig(cfg->remap, ENABLE);
	}
	
}

/**
 * @brief Function for immediately shutting off the bridge (not a braking mode) in case of faults
 * also disarms
 */
void PIOS_ESC_Off()
{
	pios_esc_dev.armed = false;
	PIOS_ESC_SetDutyCycle(0);
	PIOS_ESC_UpdateOutputs();
}

/**
 * @brief Arm the ESC.  Required for any driver outputs to go high.
 */
void PIOS_ESC_Arm()
{
	pios_esc_dev.armed = true;
}

/**
 * @brief Take a step through the state chart of modes
 */
void PIOS_ESC_NextState() 
{
	uint8_t new_state;
	
	// Simple FSM
	switch(pios_esc_dev.state) {
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
			PIOS_ESC_Off();
			return;
	}
			
	PIOS_ESC_SetState(new_state);
}

/**
 * @brief Sets the duty cycle of all PWM outputs
 */
void PIOS_ESC_SetDutyCycle(float duty_cycle)
{
	if(duty_cycle < 0 || duty_cycle > 1) {
		PIOS_ESC_Off();
		return;
	}
		
	pios_esc_dev.duty_cycle = pios_esc_dev.cfg->tim_base_init.TIM_Period * duty_cycle;
}

/**
 * @brief Sets the driver mode (which is PWM, etc)
 */
void PIOS_ESC_SetMode(enum pios_esc_mode mode)
{
	if(mode > LAST_MODE) {
		PIOS_ESC_Off();
		return;
	}
	
	pios_esc_dev.mode = mode;
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
	// This would reflect a very serious error in the logic
	if(new_state > LAST_STATE) {
		PIOS_ESC_Off();
		return;
	}
	
	pios_esc_dev.state = new_state;
	
	uint8_t first_minus;
	uint8_t first_plus;
	uint8_t second_minus;
	uint8_t second_plus;
	
	// All this logic can be made more efficient but for now making 
	// it intuitive
	
	// Based on the mode determine which gates are switching and which
	// are permanently on
	switch(pios_esc_dev.mode) {
		case ESC_MODE_LOW_ON_PWM_HIGH:
			first_minus = ESC_GATE_MODE_ON;
			first_plus = ESC_GATE_MODE_OFF;
			second_minus = ESC_GATE_MODE_OFF;
			second_plus = ESC_GATE_MODE_PWM;
			break;
		case ESC_MODE_LOW_ON_PWM_LOW:
			first_minus = ESC_GATE_MODE_ON;
			first_plus = ESC_GATE_MODE_OFF;
			second_minus = ESC_GATE_MODE_PWM_INVERT;
			second_plus = ESC_GATE_MODE_OFF;
			break;
		case ESC_MODE_LOW_ON_PWM_BOTH:
			first_minus = ESC_GATE_MODE_ON;
			first_plus = ESC_GATE_MODE_OFF;
			second_minus = ESC_GATE_MODE_PWM_INVERT;
			second_plus = ESC_GATE_MODE_PWM;
			break;
		case ESC_MODE_HIGH_ON_PWM_LOW:
			first_minus = ESC_GATE_MODE_OFF;
			first_plus = ESC_GATE_MODE_ON;
			second_minus = ESC_GATE_MODE_PWM;
			second_plus = ESC_GATE_MODE_OFF;
			break;
		case ESC_MODE_HIGH_ON_PWM_HIGH:
			first_minus = ESC_GATE_MODE_OFF;
			first_plus = ESC_GATE_MODE_ON;
			second_minus = ESC_GATE_MODE_OFF;
			second_plus = ESC_GATE_MODE_PWM_INVERT;
			break;
		case ESC_MODE_HIGH_ON_PWM_BOTH:
			first_minus = ESC_GATE_MODE_OFF;
			first_plus = ESC_GATE_MODE_ON;
			second_minus = ESC_GATE_MODE_PWM;
			second_plus = ESC_GATE_MODE_PWM_INVERT;
			break;
		default: 
			// Serious error
			PIOS_ESC_Off();
			return;
	}

	// Based on state map those gate modes to the right gates
	switch(pios_esc_dev.state) {
		case ESC_STATE_AB:
			pios_esc_dev.phase_a_minus = first_minus;
			pios_esc_dev.phase_a_plus = first_plus;			
			pios_esc_dev.phase_b_minus = second_minus;
			pios_esc_dev.phase_b_plus = second_plus;			
			pios_esc_dev.phase_c_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_c_plus = ESC_GATE_MODE_OFF;
			break;
		case ESC_STATE_AC:
			pios_esc_dev.phase_a_minus = first_minus;
			pios_esc_dev.phase_a_plus = first_plus;
			pios_esc_dev.phase_b_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_b_plus = ESC_GATE_MODE_OFF;			
			pios_esc_dev.phase_c_minus = second_minus;
			pios_esc_dev.phase_c_plus = second_plus;			
			break;
		case ESC_STATE_BA:
			pios_esc_dev.phase_a_minus = second_minus;
			pios_esc_dev.phase_a_plus = second_plus;			
			pios_esc_dev.phase_b_minus = first_minus;
			pios_esc_dev.phase_b_plus = first_plus;
			pios_esc_dev.phase_c_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_c_plus = ESC_GATE_MODE_OFF;			
			break;
		case ESC_STATE_BC:
			pios_esc_dev.phase_a_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_a_plus = ESC_GATE_MODE_OFF;			
			pios_esc_dev.phase_b_minus = first_minus;
			pios_esc_dev.phase_b_plus = first_plus;
			pios_esc_dev.phase_c_minus = second_minus;
			pios_esc_dev.phase_c_plus = second_plus;			
			break;
		case ESC_STATE_CA:
			pios_esc_dev.phase_a_minus = second_minus;
			pios_esc_dev.phase_a_plus = second_plus;			
			pios_esc_dev.phase_b_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_b_plus = ESC_GATE_MODE_OFF;			
			pios_esc_dev.phase_c_minus = first_minus;
			pios_esc_dev.phase_c_plus = first_plus;
			break;
		case ESC_STATE_CB:
			pios_esc_dev.phase_a_minus = ESC_GATE_MODE_OFF;
			pios_esc_dev.phase_a_plus = ESC_GATE_MODE_OFF;			
			pios_esc_dev.phase_b_minus = second_minus;
			pios_esc_dev.phase_b_plus = second_plus;			
			pios_esc_dev.phase_c_minus = first_minus;
			pios_esc_dev.phase_c_plus = first_plus;
			break;
		default:
			// Serious error
			PIOS_ESC_Off();
			return;			
	}
	
	PIOS_ESC_UpdateOutputs();
}

/**
 * @brief Update a single gate output
 */
static void PIOS_ESC_UpdateChannel(const struct pios_gate_channel * gate, uint8_t gate_mode) 
{
	
	uint16_t duration = pios_esc_dev.duty_cycle;
	uint8_t invert = false;
	switch(gate_mode) {
		case ESC_GATE_MODE_OFF:
			duration = 0;
			break;
		case ESC_GATE_MODE_ON:
			duration = 0xffff;
			break;
		case ESC_GATE_MODE_PWM_INVERT:
			invert = true;
		case ESC_GATE_MODE_PWM:
			duration = pios_esc_dev.duty_cycle;
			break;
		default:
			PIOS_ESC_Off();
			return;
	}	

	// If disarmed set to zero
	duration = pios_esc_dev.armed ? duration : 0;
	
	// TODO: Deal with when the output is inverted
	switch(gate->channel) {
		case TIM_Channel_1:
			TIM_SetCompare1(gate->timer, duration);
			break;
		case TIM_Channel_2:
			TIM_SetCompare2(gate->timer, duration);
			break;
		case TIM_Channel_3:
			TIM_SetCompare3(gate->timer, duration);
			break;
		case TIM_Channel_4:
			TIM_SetCompare4(gate->timer, duration);
			break;								
	}
}
/**
 * @brief Update the outputs to the gates
 */
static void PIOS_ESC_UpdateOutputs()
{
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_a_minus, pios_esc_dev.phase_a_minus);
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_a_plus, pios_esc_dev.phase_a_plus);
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_b_minus, pios_esc_dev.phase_b_minus);
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_b_plus, pios_esc_dev.phase_b_plus);
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_c_minus, pios_esc_dev.phase_c_minus);
	PIOS_ESC_UpdateChannel(&pios_esc_dev.cfg->phase_c_plus, pios_esc_dev.phase_c_plus);	
}


#endif