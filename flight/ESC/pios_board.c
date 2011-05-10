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

void PIOS_DELAY_timeout();
void PIOS_DELAY_irq_handler(void);
void TIM4_IRQHandler()
    __attribute__ ((alias("PIOS_DELAY_irq_handler")));
void PIOS_DELAY_irq_handler() {
	PIOS_DELAY_timeout();
}


#define ESC_DEFAULT_PWM_RATE 16000
#include "pios_esc_priv.h"
const struct pios_esc_cfg pios_esc_cfg = {
	.tim_base_init = {
		// Note not the same prescalar as servo
		// This is 72e6 or 24e6 / 10e6 to give a 0.1 us resolution
		.TIM_Prescaler = (PIOS_MASTER_CLOCK / 24e6) - 1,
		.TIM_ClockDivision = TIM_CKD_DIV1,
		.TIM_CounterMode = TIM_CounterMode_Up,
		.TIM_Period = ((24e6 / ESC_DEFAULT_PWM_RATE) - 1),
		.TIM_RepetitionCounter = 0x0000,
	},
	.tim_oc_init = {
		.TIM_OCMode = TIM_OCMode_PWM1,
		.TIM_OutputState = TIM_OutputState_Enable,
		.TIM_OutputNState = TIM_OutputNState_Disable,
		.TIM_Pulse = 0,		
		.TIM_OCPolarity = TIM_OCPolarity_High,
		.TIM_OCNPolarity = TIM_OCPolarity_High,
		.TIM_OCIdleState = TIM_OCIdleState_Reset,
		.TIM_OCNIdleState = TIM_OCNIdleState_Reset,
	},
	.gpio_init = {
		.GPIO_Mode = GPIO_Mode_AF_PP,
		.GPIO_Speed = GPIO_Speed_2MHz,
	},	
	.phase_a_minus = { 
		.timer = TIM2,
		.port = GPIOA,
		.channel = TIM_Channel_3,
		.pin = GPIO_Pin_2,
	},
	.phase_a_plus = {
		.timer = TIM2,
		.port = GPIOA,
		.channel = TIM_Channel_4,
		.pin = GPIO_Pin_3,
	},
	.phase_b_minus = { 
		.timer = TIM2,
		.port = GPIOA,
		.channel = TIM_Channel_1,
		.pin = GPIO_Pin_0,
	},
	.phase_b_plus = {
		.timer = TIM2,
		.port = GPIOA,
		.channel = TIM_Channel_2,
		.pin = GPIO_Pin_1,
	},
	.phase_c_minus = { /* needs to remap to alternative function */
		.timer = TIM3,
		.port = GPIOB,
		.channel = TIM_Channel_1,
		.pin = GPIO_Pin_4
	},
	.phase_c_plus = { /* needs to remap to alternative function */
		.timer = TIM3,
		.port = GPIOB,
		.channel = TIM_Channel_2,
		.pin = GPIO_Pin_5,
	},
	.remap = GPIO_PartialRemap_TIM3,
};



/*
 * ADC system
 */
#include "pios_adc_priv.h"
extern void PIOS_ADC_handler(void);
void DMA1_Channel1_IRQHandler() __attribute__ ((alias("PIOS_ADC_handler")));
// Remap the ADC DMA handler to this one
const struct pios_adc_cfg pios_adc_cfg = {
	.dma = {
		.ahb_clk  = RCC_AHBPeriph_DMA1,
		.irq = {
			.handler = PIOS_ADC_DMA_Handler,
			.flags   = (DMA1_FLAG_TC1 | DMA1_FLAG_TE1 | DMA1_FLAG_HT1 | DMA1_FLAG_GL1),
			.init    = {
				.NVIC_IRQChannel                   = DMA1_Channel1_IRQn,
				.NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGH,
				.NVIC_IRQChannelSubPriority        = 0,
				.NVIC_IRQChannelCmd                = ENABLE,
			},
		},
		.rx = {
			.channel = DMA1_Channel1,
			.init    = {
				.DMA_PeripheralBaseAddr = (uint32_t) & ADC1->DR,
				.DMA_DIR                = DMA_DIR_PeripheralSRC,
				.DMA_PeripheralInc      = DMA_PeripheralInc_Disable,
				.DMA_MemoryInc          = DMA_MemoryInc_Enable,
				.DMA_PeripheralDataSize = DMA_PeripheralDataSize_Word,
				.DMA_MemoryDataSize     = DMA_MemoryDataSize_Word,
				.DMA_Mode               = DMA_Mode_Circular,
				.DMA_Priority           = DMA_Priority_High,
				.DMA_M2M                = DMA_M2M_Disable,
			},
		}
	}, 
	.half_flag = DMA1_IT_HT1,
	.full_flag = DMA1_IT_TC1,
	.compute_downsample = false
};

struct pios_adc_dev pios_adc_devs[] = {
	{
		.cfg = &pios_adc_cfg,
		.callback_function = NULL,
	},
};

uint8_t pios_adc_num_devices = NELEMENTS(pios_adc_devs);

void PIOS_ADC_handler() {
	PIOS_ADC_DMA_Handler();
}


#if defined(PIOS_INCLUDE_USART)
#include <pios_usart_priv.h>

/*
 * DEBUG USART
 */
void PIOS_USART_debug_irq_handler(void);
void USART1_IRQHandler()
    __attribute__ ((alias("PIOS_USART_debug_irq_handler")));
const struct pios_usart_cfg pios_usart_debug_cfg = {
	.regs = USART1,
	.init = {
		 .USART_BaudRate = 230400,
		 .USART_WordLength = USART_WordLength_8b,
		 .USART_Parity = USART_Parity_No,
		 .USART_StopBits = USART_StopBits_1,
		 .USART_HardwareFlowControl = USART_HardwareFlowControl_None,
		 .USART_Mode = USART_Mode_Rx | USART_Mode_Tx,
		 },
	.irq = {
		.handler = PIOS_USART_debug_irq_handler,
		.init = {
			 .NVIC_IRQChannel = USART1_IRQn,
			 .NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_MID,
			 .NVIC_IRQChannelSubPriority = 0,
			 .NVIC_IRQChannelCmd = ENABLE,
			 },
		},
	.rx = {
	       .gpio = GPIOA,
	       .init = {
			.GPIO_Pin = GPIO_Pin_10,
			.GPIO_Speed = GPIO_Speed_2MHz,
			.GPIO_Mode = GPIO_Mode_IPU,
			},
	       },
	.tx = {
	       .gpio = GPIOA,
	       .init = {
			.GPIO_Pin = GPIO_Pin_9,
			.GPIO_Speed = GPIO_Speed_2MHz,
			.GPIO_Mode = GPIO_Mode_AF_PP,
			},
	       },
};

static uint32_t pios_usart_debug_id;
void PIOS_USART_debug_irq_handler(void)
{
	PIOS_USART_IRQ_Handler(pios_usart_debug_id);
}

#endif /* PIOS_INCLUDE_USART */

#if defined(PIOS_INCLUDE_COM)

#include <pios_com_priv.h>

#endif /* PIOS_INCLUDE_COM */

#if defined(PIOS_INCLUDE_I2C)

#include <pios_i2c_priv.h>

/*
 * I2C Adapters
 */

void PIOS_I2C_main_adapter_ev_irq_handler(void);
void PIOS_I2C_main_adapter_er_irq_handler(void);
void I2C1_EV_IRQHandler()
    __attribute__ ((alias("PIOS_I2C_main_adapter_ev_irq_handler")));
void I2C1_ER_IRQHandler()
    __attribute__ ((alias("PIOS_I2C_main_adapter_er_irq_handler")));

const struct pios_i2c_adapter_cfg pios_i2c_main_adapter_cfg = {
	.regs = I2C1,
	.init = {
		 .I2C_Mode = I2C_Mode_I2C,
		 .I2C_OwnAddress1 = 0,
		 .I2C_Ack = I2C_Ack_Enable,
		 .I2C_AcknowledgedAddress = I2C_AcknowledgedAddress_7bit,
		 .I2C_DutyCycle = I2C_DutyCycle_2,
		 .I2C_ClockSpeed = 200000,	/* bits/s */
		 },
	.transfer_timeout_ms = 50,
	.scl = {
		.gpio = GPIOB,
		.init = {
			 .GPIO_Pin = GPIO_Pin_6,
			 .GPIO_Speed = GPIO_Speed_10MHz,
			 .GPIO_Mode = GPIO_Mode_AF_OD,
			 },
		},
	.sda = {
		.gpio = GPIOB,
		.init = {
			 .GPIO_Pin = GPIO_Pin_7,
			 .GPIO_Speed = GPIO_Speed_10MHz,
			 .GPIO_Mode = GPIO_Mode_AF_OD,
			 },
		},
	.event = {
		  .handler = PIOS_I2C_main_adapter_ev_irq_handler,
		  .flags = 0,	/* FIXME: check this */
		  .init = {
			   .NVIC_IRQChannel = I2C1_EV_IRQn,
			   .NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGHEST,
			   .NVIC_IRQChannelSubPriority = 0,
			   .NVIC_IRQChannelCmd = ENABLE,
			   },
		  },
	.error = {
		  .handler = PIOS_I2C_main_adapter_er_irq_handler,
		  .flags = 0,	/* FIXME: check this */
		  .init = {
			   .NVIC_IRQChannel = I2C1_ER_IRQn,
			   .NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGHEST,
			   .NVIC_IRQChannelSubPriority = 0,
			   .NVIC_IRQChannelCmd = ENABLE,
			   },
		  },
};

uint32_t pios_i2c_main_adapter_id;
void PIOS_I2C_main_adapter_ev_irq_handler(void)
{
	/* Call into the generic code to handle the IRQ for this specific device */
	PIOS_I2C_EV_IRQ_Handler(pios_i2c_main_adapter_id);
}

void PIOS_I2C_main_adapter_er_irq_handler(void)
{
	/* Call into the generic code to handle the IRQ for this specific device */
	PIOS_I2C_ER_IRQ_Handler(pios_i2c_main_adapter_id);
}

#endif /* PIOS_INCLUDE_I2C */


extern const struct pios_com_driver pios_usart_com_driver;

uint32_t pios_com_debug_id;
uint8_t adc_fifo_buf[sizeof(float) * 6 * 4] __attribute__ ((aligned(4)));    // align to 32-bit to try and provide speed improvement

/**
 * PIOS_Board_Init()
 * initializes all the core subsystems on this specific hardware
 * called from System/openpilot.c
 */
void PIOS_Board_Init(void) {
	/* Brings up System using CMSIS functions, enables the LEDs. */
	PIOS_SYS_Init();

	PIOS_LED_Off(LED1);
	PIOS_LED_Off(LED2);
	PIOS_LED_Off(LED3);

	/* Delay system */
	PIOS_DELAY_Init();
	
	PIOS_GPIO_Init();
	
	/* Bring up ADC for sensing BEMF */
	PIOS_ADC_Init();
	PIOS_ADC_Config(5);
	
	/* Remap AFIO pin */
	GPIO_PinRemapConfig( GPIO_Remap_SWJ_NoJTRST, ENABLE);
	PIOS_ESC_Init(&pios_esc_cfg);

	/* Communication system */
	if (PIOS_USART_Init(&pios_usart_debug_id, &pios_usart_debug_cfg)) {
		PIOS_DEBUG_Assert(0);
	}
	if (PIOS_COM_Init(&pios_com_debug_id, &pios_usart_com_driver, pios_usart_debug_id)) {
		PIOS_DEBUG_Assert(0);
	}
}

