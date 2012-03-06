/**
  ******************************************************************************
  * @file swuart.h
  * @brief This file contains all user and system definition of uart software.
  * @author STMicroelectronics - MCD Application Team
  * @version V1.0.0
  * @date 10/13/2008
  ******************************************************************************
  *
  * THE PRESENT SOFTWARE WHICH IS FOR GUIDANCE ONLY AIMS AT PROVIDING CUSTOMERS
  * WITH CODING INFORMATION REGARDING THEIR PRODUCTS IN ORDER FOR THEM TO SAVE
  * TIME. AS A RESULT, STMICROELECTRONICS SHALL NOT BE HELD LIABLE FOR ANY
  * DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES WITH RESPECT TO ANY CLAIMS ARISING
  * FROM THE CONTENT OF SUCH SOFTWARE AND/OR THE USE MADE BY CUSTOMERS OF THE
  * CODING INFORMATION CONTAINED HEREIN IN CONNECTION WITH THEIR PRODUCTS.
  *
  * <h2><center>&copy; COPYRIGHT 2008 STMicroelectronics</center></h2>
  * @image html logo.bmp
  ******************************************************************************
  */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __SW_UART_H
#define __SW_UART_H

// ---------------------------------------------------------
// --------------- user macros definition ------------------
// ---------------------------------------------------------
#ifndef _COSMIC_
#define _Bool   bit
#endif

#define SWUART_TRANSMIT_USED
#define SWUART_RECEIVE_USED
//#define TEST_PIN_USED

#define	DATA_LENGTH		8		// number of data bits in protocol
#define	STOP_BITS		1		// number of stop bits in protocol
//#define  PARITY
// if PARITY flag si defined the most significant data
// 			bit is considered to be parity bit

// user HW definition -------------------------------------
// Macros in following section should be carefully completed by user 
// in according with hardware connection of Tx and Rx pins, used
// protocol and speed (Rx pin is supposed with capture capability,
// Tx pin is any with output capability)
//----------------------------------------------------------

#ifdef SWUART_TRANSMIT_USED
#define  UART_TxPORT	(GPIOD)
#define	UART_TxPIN	(GPIO_PIN_2)
#endif

#ifdef SWUART_RECEIVE_USED
#define  UART_RxPORT	(GPIOD)
#define	UART_RxPIN	(GPIO_PIN_0)
#endif

#ifdef TEST_PIN_USED
#define  UART_TestPORT	(GPIOD)
#define	UART_TestPIN	(GPIO_PIN_4)
#endif

#define	OCInit_TypeDef		TIM3_OCInit_TypeDef TIM3_OCInitStruct

/* uncomment the next line to switch the system clock from HSI (16 Mhz) to HSE (24 Mhz) */
#define TO_HSE

//	Here user must ensure calling service interrupt routine at exact intervals 
//	of 1/2 bit of the dedicated baud rate (see also TO_HSE switch above)
// TIM3 overflow period must be set in dependance on fmaster clock selected:
//                                     fmaster= fcpu clock:   16 / 24 Mhz
// 	for 9600Bd speed set overflow every 52.08us - Period=  833 / 1250
//   for 19200Bd speed set overflow every 26,04us - Period=  417 / 625
//   for 57600Bd speed set overflow every  8,68us - Period=  139*/ 208
// *) for half duplex only!
// next macro initializes and enable TIM3 base and interrupt system, divider
// register /1, auto reload register, enable interrupt on update, update on overflow
#define init_UART_timing	{\
	TIM3_TimeBaseInit(TIM3_PRESCALER_1, 208);\
	TIM3_ITConfig(TIM3_IT_UPDATE, ENABLE);\
	TIM3_UpdateRequestConfig(TIM3_UPDATESOURCE_REGULAR);\
	TIM3_Cmd(ENABLE);\
}

// input capture system initialization
//filter (4 samples), IC2 mapped on TI2FP2, capture on faling edge
#define	init_ic_setting	{\
	TIM3->CCER1&= ~TIM3_CCER1_CC2E;\
	TIM3->CCMR2= (((2<<4) & TIM3_CCMR_ICxF) | ((1<<0) & TIM3_CCMR_CCxS));\
	TIM3->CCER1|= TIM3_CCER1_CC2P | TIM3_CCER1_CC2E;\
}
// output compare system initialization
// frozen (timing only) mode, preload regs disable, no output
#define  init_oc_setting {\
	TIM3->CCER1&= ~TIM3_CCER1_CC2E;\
	TIM3->CCMR2= 0;\
}

#define	clear_owerflow_flag	{ TIM3->SR1 &= ~TIM3_SR1_UIF; }
#define	clear_cc_flag			{ TIM3->SR1 &= ~TIM3_SR1_CC2IF;	}
#define  enable_cc_interrupt	{ TIM3->IER |=  TIM3_IER_CC2IE; }
#define  disable_cc_interrupt	{ TIM3->IER &= ~TIM3_IER_CC2IE; }

#define	enable_IC_system		{ /* clear_cc_flag;*/ init_ic_setting; enable_cc_interrupt; }
#define	disable_IC_system	 {\
	disable_cc_interrupt;\
	TIM3->CCER1 &= (u8)(~TIM3_CCER1_CC2E);\
	TIM3->CCMR2&=~ TIM3_CCMR_CCxS;\
};

#define	enable_OC_system		{ /* clear_cc_flag;*/ init_oc_setting; enable_cc_interrupt; }
#define	disable_OC_system	  	  disable_cc_interrupt

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

// bit manipulation definition
#ifdef SWUART_RECEIVE_USED
#define Rx_test	(UART_RxPORT->IDR & UART_RxPIN)
#endif
#ifdef SWUART_TRANSMIT_USED
#define set_Tx		(UART_TxPORT->ODR |= UART_TxPIN)
#define clr_Tx		(UART_TxPORT->ODR &=~UART_TxPIN)
#endif
#ifdef TEST_PIN_USED
#define set_Test_Pin		(UART_TestPORT->ODR |= UART_TestPIN)
#define clr_Test_Pin		(UART_TestPORT->ODR &=~UART_TestPIN)
#endif

// initial definition of HW pins Tx as output push-pull, Rx as input floating
 
#ifdef SWUART_TRANSMIT_USED
#define init_UART_Tx_pin	{ GPIO_Init(UART_TxPORT, UART_TxPIN, GPIO_MODE_OUT_PP_HIGH_SLOW); }
#endif
#ifdef SWUART_RECEIVE_USED
#define init_UART_Rx_pin	{ GPIO_Init(UART_RxPORT, UART_RxPIN, GPIO_MODE_IN_FL_NO_IT); }
#endif
#ifdef TEST_PIN_USED
#define init_UART_Test_pin	{ GPIO_Init(UART_TestPORT, UART_TestPIN, GPIO_MODE_OUT_PP_LOW_FAST); }
#endif

/* Exported constants --------------------------------------------------------*/
// status masks
#define	transmit_in_progress			0x80
#define	transmit_data_reg_empty		0x40
#define	receive_9th_data_bit			0x20
#define	receive_in_progress			0x10

#define	receive_buffer_overflow		0x08 // low nibble corresponds to error return value
#define	receive_frame_error			0x04
#define	receive_noise_error			0x02
#define	receive_buffer_full			0x01

// error return codes
#define OVFL_ERROR	8
#define FRAME_ERROR	4
#define NOISE_ERROR	2
#define RX_BUFF_FULL	1

        
/* Exported macros ------------------------------------------------------------*/
#define test_status(a)		(UART_sts & a)
#define set_status(a)		(UART_sts |= a)
#define clr_status(a)		(UART_sts &=~a)

#define uart_receive_enable 	{ Rx_bit= Rx_phase= 0; disable_OC_system; clear_cc_flag;enable_IC_system; }
#define uart_receive_disable 	{ disable_IC_system; disable_OC_system; clr_status(receive_in_progress); }

/* Exported variables --------------------------------------------------------*/
extern _Bool Rx_phase;
extern _Bool Tx_phase;
#ifdef PARITY
extern _Bool Rx_parity;
extern _Bool Tx_parity;
#else
#ifdef BIT9
extern _Bool Rx_bit9;
extern _Bool Tx_bit9;
#endif
#endif
extern u8 Rx_bit,
			 Rx_samp,
			 Tx_bit,
			 Rx_buff,
			 Rx_data,
			 Tx_data,
			 UART_sts;
			 
/* Exported functions ------------------------------------------------------- */
void uart_init(void);
u8 uart_send(u8 b);
u8 uart_read(u8 *b);
void uart_Tx_timing(void);
void uart_Rx_timing(void);
#endif

/******************* (C) COPYRIGHT 2008 STMicroelectronics *****END OF FILE****/
