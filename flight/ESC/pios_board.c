/**
 ******************************************************************************
 *
 * @file       pios_board.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      Defines board specific static initializers for hardware for the AHRS board.
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

#include <pios.h>

#include "board_hw_defs.c"

/**
 * PIOS_Board_Init()
 * initializes all the core subsystems on this specific hardware
 * called from System/openpilot.c
 */
void PIOS_Board_Init(void) {
	/* Brings up System using CMSIS functions, enables the LEDs. */
	PIOS_SYS_Init();

	/* Delay system */
	PIOS_DELAY_Init();
	
	PIOS_LED_Init(&pios_led_cfg);
	
	PIOS_GPIO_Init();

	/* Bring up ADC for sensing BEMF */
	PIOS_ADC_Init();
	PIOS_ADC_Config(1);

	/* Remap AFIO pin */
//	GPIO_PinRemapConfig( GPIO_Remap_SWJ_NoJTRST, ENABLE);
	PIOS_ESC_Init(&pios_esc_cfg);

	PIOS_TIM_InitClock(&tim_4_cfg);
#if defined(PIOS_INCLUDE_PWM)
	uint32_t pios_pwm_id;
	PIOS_PWM_Init(&pios_pwm_id, &pios_pwm_cfg);
	
	uint32_t pios_pwm_rcvr_id;
	if (PIOS_RCVR_Init(&pios_pwm_rcvr_id, &pios_pwm_rcvr_driver, pios_pwm_id)) {
		PIOS_Assert(0);
	}
	pios_rcvr_group_map[0] = pios_pwm_rcvr_id;
#endif

	/* Communication system */
	if (PIOS_USART_Init(&pios_usart_debug_id, &pios_usart_debug_cfg)) {
		PIOS_DEBUG_Assert(0);
	}

	if (PIOS_COM_Init(&pios_com_debug_id, &pios_usart_com_driver, pios_usart_debug_id,
					  uart_debug_rx_buffer,UART_DEBUG_RX_LEN,
					  uart_debug_tx_buffer,UART_DEBUG_TX_LEN)) {
		PIOS_DEBUG_Assert(0);
	}
	
	PIOS_WDG_Init();
}

