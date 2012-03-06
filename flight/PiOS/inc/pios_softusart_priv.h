/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup   PIOS_SOFTUSART USART Functions
 * @brief PIOS interface for emulated USART port
 * @{
 *
 * @file       pios_softusart_priv.h
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2012.
 * @brief      USART private definitions.
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

// Adopted from ST example AN2781

#ifndef PIOS_SOFTUSART_H
#define PIOS_SOFTUSART_H

extern const struct pios_com_driver pios_softusart_com_driver;

struct pios_softusart_cfg {
	TIM_ICInitTypeDef tim_ic_init;
	struct pios_tim_channel rx;
	struct pios_tim_channel tx;
	bool half_duplex;
};

extern int32_t PIOS_SOFTUSART_Init(uint32_t * softusart_id, const struct pios_softusart_cfg * cfg);
extern const struct pios_softusart_cfg * PIOS_SOFTUSART_GetConfig(uint32_t softusart_id);

#define	DATA_LENGTH		8		// number of data bits in protocol
#define	STOP_BITS		1		// number of stop bits in protocol
//#define  PARITY
// if PARITY flag si defined the most significant data
// 			bit is considered to be parity bit

// ---------------------------------------------------------
// ------------- private macros definition -----------------
// ---------------------------------------------------------

#if (DATA_LENGTH == 9)
#ifndef PARITY
#define BIT9					// definition of 9-th data bit flag
#define	set_Tx_bit9		{ Tx_bit9= 1; }
#define	clr_Tx_bit9		{ Tx_bit9= 0; }
#endif
#endif

#if (DATA_LENGTH > 9)				// checking the parameters range
#error DATA LENGTH SHOULD NOT BE LONGER THAN NINE BITS
#endif
#if (DATA_LENGTH < 1)
#error DATA LENGTH SHOULD NOT BE SHORTER THAN ONE BIT
#endif
#if (STOP_BITS > 2)
#error NO MORE THAN TWO STOP BITS SHOULD BE DEFINED
#endif
#if (STOP_BITS < 1)
#error AT LEAST ONE STOP BIT SHOULD BE DEFINED
#endif

/* Exported constants --------------------------------------------------------*/
// status masks
#define	TRANSMIT_IN_PROGRESS        0x80
#define	TRANSMIT_DATA_REG_EMPTY     0x40
#define	RECEIVED_9TH_DATA_BIT       0x20
#define	RECEIVE_IN_PROGRESS         0x10
#define	RECEIVE_BUFFER_OVERFLOW     0x08 // low nibble corresponds to error return value
#define	RECEIVE_FRAME_ERROR         0x04
#define	RECEIVE_NOISE_ERROR         0x02


// error return codes
#define OVFL_ERROR	8
#define FRAME_ERROR	4
#define NOISE_ERROR	2
#define RX_BUFF_FULL	1

//#define uart_receive_enable 	{ Rx_bit= Rx_phase= 0; disable_OC_system; clear_cc_flag;enable_IC_system; }
//#define uart_receive_disable 	{ disable_IC_system; disable_OC_system; clr_status(receive_in_progress); }

/* Exported variables --------------------------------------------------------*/
#endif
