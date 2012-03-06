/**
 ******************************************************************************
 * @addtogroup PIOS PIOS Core hardware abstraction layer
 * @{
 * @addtogroup   PIOS_PWM PWM Input Functions
 * @brief		Code to measure with PWM input
 * @{
 *
 * @file       pios_pwm.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      PWM Input functions (STM32 dependent)
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

#if defined(PIOS_INCLUDE_SOFTUSART)

#include "pios_softusart_priv.h"

/* Provide a COM driver */
static void PIOS_SOFTUSART_ChangeBaud(uint32_t usart_id, uint32_t baud);
static void PIOS_SOFTUSART_RegisterRxCallback(uint32_t usart_id, pios_com_callback rx_in_cb, uint32_t context);
static void PIOS_SOFTUSART_RegisterTxCallback(uint32_t usart_id, pios_com_callback tx_out_cb, uint32_t context);
static void PIOS_SOFTUSART_TxStart(uint32_t usart_id, uint16_t tx_bytes_avail);
static void PIOS_SOFTUSART_RxStart(uint32_t usart_id, uint16_t rx_bytes_avail);

const struct pios_com_driver pios_usart_com_driver = {
	.set_baud   = PIOS_SOFTUSART_ChangeBaud,
	.tx_start   = PIOS_SOFTUSART_RegisterRxCallback,
	.rx_start   = PIOS_SOFTUSART_RegisterTxCallback,
	.bind_tx_cb = PIOS_SOFTUSART_TxStart,
	.bind_rx_cb = PIOS_SOFTUSART_RxStart,
};

/* Local Variables */
enum pios_softusart_dev_magic {
	PIOS_SOFTUSART_DEV_MAGIC = 0xab30293c,
};

const uint8_t MSK_TAB[9]= { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0};

struct pios_softusart_dev {
	enum pios_softusart_dev_magic     magic;
	const struct pios_softusart_cfg  *cfg;
	
	// Comm variables
	pios_com_callback rx_in_cb;
	uint32_t rx_in_context;
	pios_com_callback tx_out_cb;
	uint32_t tx_out_context;

	// Communication variables
	bool rx_phase;     // phase of received bit [0-1] (edge, middle)
	bool tx_phase;     // phase of transmited bit [0-1] (edge, middle)
	bool rx_parity;    // received parity [0-1]
	bool tx_parity;    // transmited parity [0-1]
	bool rx_pin9;      // received 9-th data bit [0-1]
	bool tx_pin9;      // transmited 9-th data bit [0-1]
	uint8_t rx_bit;    // counter of received bits [0-11]
	uint8_t tx_bit;    // counter of transmited bits [0-11]
	uint8_t rx_samp;   // register of samples [0-3]
	uint8_t rx_buff;
	uint8_t rx_data;   // received byte register
	uint8_t tx_data;   // transmited byte register
	uint8_t status;	   // UART status register (1= active state)
};


static bool PIOS_SOFTUSART_TestStatus(struct pios_softusart_dev *dev, flag);
static void PIOS_SOFTUSART_SetStatus(struct pios_softusart_dev *dev, flag);
static void PIOS_SOFTUSART_SetStatus(struct pios_softusart_dev *dev, flag);

static bool PIOS_SOFTUSART_validate(struct pios_softusart_dev * pwm_dev)
{
	return (pwm_dev->magic == PIOS_SOFTUSART_DEV_MAGIC);
}

#if defined(PIOS_INCLUDE_FREERTOS)
static struct pios_softusart_dev * PIOS_SOFTUSART_alloc(void)
{
	struct pios_softusart_dev *pwm_softusart;

	pwm_softusart = (struct pios_softusart_dev *)pvPortMalloc(sizeof(*pwm_softusart));
	if (!pwm_softusart) return(NULL);

	pwm_softusart->magic = PIOS_SOFTUSART_DEV_MAGIC;
	return(pwm_softusart);
}
#else
static struct pios_softusart_dev pios_softusart_devs[PIOS_SOFTUSART_MAX_DEVS];
static uint8_t pios_softusart_num_devs;
static struct pios_softusart_dev * PIOS_SOFTUSART_alloc(void)
{
	struct pios_softusart_dev *pwm_softusart;

	if (pios_softusart_num_devs >= PIOS_SOFTUSART_MAX_DEVS) {
		return (NULL);
	}

	pwm_softusart = &pios_pwm_devs[pios_softusart_num_devs++];
	pwm_softusart->magic = PIOS_SOFTUSART_DEV_MAGIC;

	return (pwm_softusart);
}
#endif

static void PIOS_SOFTUSART_tim_overflow_cb (uint32_t id, uint32_t context, uint8_t channel, uint16_t count);
static void PIOS_SOFTUSART_tim_edge_cb (uint32_t id, uint32_t context, uint8_t channel, uint16_t count);
const static struct pios_tim_callbacks tim_callbacks = {
	.overflow = PIOS_SOFTUSART_tim_overflow_cb,
	.edge     = PIOS_SOFTUSART_tim_edge_cb,
};

/**
* Initialises all the pins
*/
int32_t PIOS_SOFTUSART_Init(uint32_t *softusart_id, const struct pios_softusart_cfg *cfg)
{
	PIOS_DEBUG_Assert(pwm_id);
	PIOS_DEBUG_Assert(cfg);

	struct pios_softusart_dev *pwm_softusart;

	pwm_softusart = (struct pios_softusart_dev *) PIOS_SOFTUSART_alloc();
	if (!pwm_softusart) goto out_fail;

	/* Bind the configuration to the device instance */
	pwm_softusart->cfg = cfg;

	uint32_t tim_id;
	if (PIOS_TIM_InitChannels(&tim_id, cfg->channels, cfg->num_channels, &tim_callbacks, (uint32_t)pwm_softusart)) {
		return -1;
	}

	/* Configure the channels to be in capture/compare mode */
	for (uint8_t i = 0; i < cfg->num_channels; i++) {
		const struct pios_tim_channel * chan = &cfg->channels[i];

		/* Configure timer for input capture */
		TIM_ICInitTypeDef TIM_ICInitStructure = cfg->tim_ic_init;
		TIM_ICInitStructure.TIM_Channel = chan->timer_chan;
		TIM_ICInit(chan->timer, &TIM_ICInitStructure);
		
		/* Enable the Capture Compare Interrupt Request */
		switch (chan->timer_chan) {
		case TIM_Channel_1:
			TIM_ITConfig(chan->timer, TIM_IT_CC1, ENABLE);
			break;
		case TIM_Channel_2:
			TIM_ITConfig(chan->timer, TIM_IT_CC2, ENABLE);
			break;
		case TIM_Channel_3:
			TIM_ITConfig(chan->timer, TIM_IT_CC3, ENABLE);
			break;
		case TIM_Channel_4:
			TIM_ITConfig(chan->timer, TIM_IT_CC4, ENABLE);
			break;
		}

		// Need the update event for that timer to detect timeouts
		TIM_ITConfig(chan->timer, TIM_IT_Update, ENABLE);
	}

	*softusart_id = (uint32_t) pwm_softusart;

	return (0);

out_fail:
	return (-1);
}

/**
 * @brief Check a status flag
 */
static bool PIOS_SOFTUSART_TestStatus(struct pios_softusart_dev *dev, flag)
{
	return dev->status & flag;
}

/**
 * @brief Set a status flag
 */
static void PIOS_SOFTUSART_SetStatus(struct pios_softusart_dev *dev, flag)
{
	dev->status |= flag;
}

/**
 * @brief Clear a status flag
 */
static void PIOS_SOFTUSART_ClrStatus(struct pios_softusart_dev *dev, flag)
{
	dev->status &= ~flag;
}

/**
 * @brief Set the baud rate
 */
static void PIOS_SOFTUSART_ChangeBaud(uint32_t usart_id, uint32_t baud)
{

}

/**
 * @brief Set the callback into the general com driver when a byte is received
 */
static void PIOS_SOFTUSART_RegisterRxCallback(uint32_t usart_id, pios_com_callback rx_in_cb, uint32_t context)
{
	struct pios_softusart_dev *softusart_dev = (struct pios_softusart_dev *) usart_id;
	bool valid = PIOS_SOFTUSART_validate(softusart_dev);
	PIOS_Assert(valid);
	
	/* 
	 * Order is important in these assignments since ISR uses _cb
	 * field to determine if it's ok to dereference _cb and _context
	 */
	softusart_dev->rx_in_context = context;
	softusart_dev->rx_in_cb = rx_in_cb;	
}

/**
 * @brief Set the callback into the general com driver when a byte should be transmitted
 */
static void PIOS_SOFTUSART_RegisterTxCallback(uint32_t usart_id, pios_com_callback tx_out_cb, uint32_t context)
{
	struct pios_softusart_dev *softusart_dev = (struct pios_softusart_dev *) usart_id;
	bool valid = PIOS_SOFTUSART_validate(softusart_dev);
	PIOS_Assert(valid);
	
	/* 
	 * Order is important in these assignments since ISR uses _cb
	 * field to determine if it's ok to dereference _cb and _context
	 */
	softusart_dev->tx_out_context = context;
	softusart_dev->tx_out_cb = tx_out_cb;
}

/**
 * @brief Start transmission
 * @param[in] usart_id the context for the usart device
 * @parma[in] tx_bytes_avail how many bytes are available to transmit
 */
static void PIOS_SOFTUSART_TxStart(uint32_t usart_id, uint16_t tx_bytes_avail)
{
	struct pios_softusart_dev *softusart_dev = (struct pios_softusart_dev *) usart_id;

	if(!PIOS_SOFTUSART_TestStatus(softusart_dev, TRANSMIT_IN_PROGRESS) ||
	    PIOS_SOFTUSART_TestStatus(softusart_dev, TRANSMIT_DATA_REG_EMPTY)) {
		//YES - initiate sending procedure
		softusart_dev->tx_data = b;             
		PIOS_SOFTUSART_ClrStatus(softusart_dev, TRANSMIT_DATA_REG_EMPTY);
		if(!PIOS_SOFTUSART_TestStatus(softusart_dev, TRANSMIT_IN_PROGRESS)) {
			softusart_dev->tx_phase = 0;
			softusart_dev->tx_bit = 0;
			PIOS_SOFTUSART_SetStatus(softusart_dev, TRANSMIT_IN_PROGRESS);
		};
		return true;
	}
	else
		return false;  //NO - no action (transmition in progress)
}

static void PIOS_SOFTUSART_RxStart(uint32_t usart_id, uint16_t rx_bytes_avail)
{
}

void uart_Tx_timing(void);
void uart_Rx_timing(void);

/**
 * @brief When this occurs determine whether to set output high or low
 */
static void PIOS_SOFTUSART_tim_overflow_cb (uint32_t tim_id, uint32_t context, uint8_t channel, uint16_t count)
{
	struct pios_softusart_dev *softusart_dev = (struct pios_softusart_dev *) context;

	if (!PIOS_SOFTUSART_validate(softusart_dev)) {
		/* Invalid device specified */
		return;
	}

	if (channel >= softusart_dev->cfg->num_channels) {
		/* Channel out of range */
		return;
	}

	clear_owerflow_flag;
	if(softusart_dev->tx_phase) {
		if(PIOS_SOFTUSART_TestStatus(softusart_dev, TRANSMIT_IN_PROGRESS)) { // edge of current bit (no service for middle)
			switch(softusart_dev->tx_bit) {     // begin of bit transmition
				case 0:
					softusart_dev->tx_bit9 = 0; //start bit transmition
#ifdef PARITY
					softusart_dev->tx_parity = 0;
#endif
					break;
#ifdef PARITY
				case DATA_LENGTH:		
					if(softusart->tx_parity) 
						set_Tx;
					else
						clr_Tx;
					break;
#else
#ifdef BIT9
				case DATA_LENGTH:		
					if(softusart_dev->tx_bit9)
						set_Tx;
					else
						clr_Tx;
					break;
#endif
#endif
				case DATA_LENGTH+1:	
					PIOS_SOFTUSART_SetStatus(softusart_dev,TRANSMIT_DATA_REG_EMPTY);
				case DATA_LENGTH+2:	
					set_Tx;	//stop bit(s) transmition
					break;
#ifdef PARITY
				default:	
					if(softusart_dev->tx_data & MSK_TAB[softusart_dev->tx_bit-1])
					{ 
						set_Tx; 
						softusart_dev->tx_parity = ~softusart_dev->tx_parity;
					}
					else
						clr_Tx; // parity transmition
#else
				default:	
					if(softusart_dev->tx_data & MSK_TAB[softusart_dev->tx_bit-1])
						set_Tx; //  data bits
					else
						clr_Tx; // transmition
#endif
			};
			if(softusart_dev->tx_bit >= DATA_LENGTH + STOP_BITS) {
			    // end of current transmited bit
				if(PIOS_SOFTUSART_TestStatus(softusart_dev, TRANSMIT_DATA_REG_EMPTY)) {	// end of transmition
					PIOS_SOFTUSART_ClrStatus(softusart_dev, TRANSMIT_IN_PROGRESS);      // no next Tx request
				}
				else
					softusart_dev->tx_bit= 0;     // next byte is buffered - continue
			}                                     // with transmition
			else
				++softusart_dev->tx_bit;          // next bit to transmit			
		};
	};
	softusart_dev->tx_phase = ~softusart_dev->tx_phase;
	return;
}

static void PIOS_SOFTUSART_tim_edge_cb (uint32_t tim_id, uint32_t context, uint8_t chan_idx, uint16_t count)
{
	/* Recover our device context */
		struct pios_softusart_dev *softusart_dev = (struct pios_softusart_dev *) context;

	if (!PIOS_SOFTUSART_validate(softusart_dev)) {
		/* Invalid device specified */
		return;
	}

	if (channel >= softusart_dev->cfg->num_channels) {
		/* Channel out of range */
		return;
	}

	clear_cc_flag;
	if(PIOS_SOFTUSART_TestStatus(softusart_dev, RECEIVE_IN_PROGRESS)) {
		if(!softusart_dev->rx_phase) {      // receive is in progres now
			softusart_dev->rx_samp= 0;      // middle of current bit, checking samples

			if(softusart_dev->rx_test)	++softusart_dev->rx_samp; // sampling in the middle of current bit
			if(softusart_dev->rx_test)	++softusart_dev->rx_samp;
			if(softusart_dev->rx_test)	++softusart_dev->rx_samp;
			
			if(softusart_dev->rx_bit == 0) {
				if(softusart_dev->rx_samp == 0) {  // start bit!
					softusart_dev->rx_bit = 1;      // correctly received, continue
					softusart_dev->rx_buff = 0;
				} else {					// noise in start bit, find next one
					disable_OC_system;
					enable_IC_system;
					PIOS_SOFTUSART_ClrStatus(softusart_dev, RECEIVE_IN_PROGRESS);
				};
			}
			else {
				switch(softusart_dev->rx_samp) {		// any other bit, results?
					case 1:	
						PIOS_SOFTUSART_SetStatus(softusart_dev, RECEIVE_NOISE_ERROR);
						break;	// noise in middle samples, "0" received
					case 2: 	
						PIOS_SOFTUSART_SetStatus(softusart_dev, RECEIVE_NOISE_ERROR);
						// noise in middle samples, "1" received
#ifdef PARITY
					case 3:	
						if(softusart_dev->rx_bit < DATA_LENGTH)
							softusart_dev->rx_buff|= MSK_TAB[softusart_dev->rx_bit-1];
						if(softusart_dev->rx_bit <= DATA_LENGTH)
							softusart_dev->rx_parity = ~softusart_dev->rx_parity;
#else
#ifdef BIT9
					case 3:	
						if(softusart_dev->rx_bit < DATA_LENGTH)
							softusart_dev->rx_buff |= MSK_TAB[softusart_dev->rx_bit-1];
						if(softusart_dev->rx_bit == DATA_LENGTH)
							softusart_dev->rx_bit9=  1;
#else
					case 3:
						if(softusart_dev->rx_bit <= DATA_LENGTH)
							softusart_dev->rx_buff |= MSK_TAB[softusart_dev->rx_bit-1];
#endif
#endif
						break;	// "1" correctly received
				};
				if(softusart_dev->rx_bit > DATA_LENGTH) {
#ifdef PARITY
					// stop bit(s) are received, results?
					if(softusart_dev->rx_samp != 3  || softusart_dev->rx_parity)
#else
						if(softusart_dev->rx_samp != 3)	// stop bit(s) are received, results?
#endif
							PIOS_SOFTUSART_SetStatus(softusart_dev,RECEIVE_FRAME_ERROR);// noise in stop bit or parity error
					if(softusart_dev->rx_bit >= DATA_LENGTH + STOP_BITS) {
						if(!PIOS_SOFTUSART_TestStatus(softusart_dev,RECEIVE_BUFFER_FULL)) { // end of receive
							softusart_dev->rx_data = softusart_dev->rx_buff; // new byte in buffer
#ifdef BIT9
							if(softusart_dev->rx_bit9)
								PIOS_SOFTUSART_SetStatus(softusart_dev,RECEIVED_9TH_DATA_BIT);
							else
								PIOS_SOFTUSART_ClrStatus(softusart_dev,RECEIVED_9TH_DATA_BIT);
#endif
							PIOS_SOFTUSART_SetStatus(softusart_dev,RECEIVE_BUFFER_FULL);		 
						}
						else
							PIOS_SOFTUSART_SetStatus(softusart_dev,RECEIVE_BUFFER_OVERFLOW); // data overflow!
						
						disable_OC_system;		// init next byte receive
						PIOS_SOFTUSART_ClrStatus(softusart_dev,RECEIVE_IN_PROGRESS);
						enable_IC_system;
					}
					else
						++softusart_dev->rx_bit;
				}
				else
					++softusart_dev->rx_bit;			// init next data bit receive
			}
		}
		softusart_dev->rx_phase = ~softusart_dev->rx_phase;
	}
	else {								// receive is not in progres yet
		disable_IC_system;			    // IC interrupt - begin of start bit detected
		enable_OC_system;				//	OC interrupt period as a conseqence of start bit falling edge
		PIOS_SOFTUSART_SetStatus(softusart_dev,RECEIVE_IN_PROGRESS);	// receive byte initialization
		softusart_dev->rx_bit = 0;
		softusart_dev->rx_phase = 0;
#ifdef PARITY
		softusart_dev->rx_parity = 0;
#else
#ifdef BIT9
		softusart_dev->rx_bit9 = 0;
#endif
#endif
	};
}

/**
 * @brief Sending byte initialization if it could be buffered
 * @par Parameters:
 * u8 byte to be transmitted
 * @retval TRUE (byte is sent) or FALSE (byte could not be buffered)
 */
uint8_t uart_send(uint8_t b) {	//transmition/bufferring possible?
}
//----------------------------------------------------------
/**
 * @brief Timing of uart Tx - all bits transmit proccess control
 * @par Parameters:
 * None
 * @retval None
 */
void uart_Tx_timing(void) { 
	clear_owerflow_flag;
	if(Tx_phase) {
		if(test_status(transmit_in_progress)) {	// edge of current bit (no service for middle)
			switch(Tx_bit) {			// begin of bit transmition
				case 0:					clr_Tx;	//start bit transmition
#ifdef PARITY
					Tx_parity= 0;
#endif
					break;
#ifdef PARITY
				case DATA_LENGTH:		if(Tx_parity) set_Tx;
				else          clr_Tx;
					break;
#else
#ifdef BIT9
				case DATA_LENGTH:		if(Tx_bit9) set_Tx;
				else        clr_Tx;
					break;
#endif
#endif
				case DATA_LENGTH+1:	set_status(transmit_data_reg_empty);
				case DATA_LENGTH+2:	set_Tx;	//stop bit(s) transmition
					break;
#ifdef PARITY
				default:	if(Tx_data & MSK_TAB[Tx_bit-1]) { set_Tx; Tx_parity= ~Tx_parity; }
				else										 clr_Tx; // parity transmition
#else
				default:	if(Tx_data & MSK_TAB[Tx_bit-1])	set_Tx; //  data bits
				else										clr_Tx; // transmition
#endif
			};
			if(Tx_bit >= DATA_LENGTH + STOP_BITS) {	// end of current transmited bit
				if(test_status(transmit_data_reg_empty)) {	// end of transmition
					clr_status(transmit_in_progress);				// no next Tx request
				}
				else
					Tx_bit= 0;							// next byte is buffered - continue
			}												//			with transmition
			else
				++Tx_bit;				// next bit to transmit			
		};
	};
	Tx_phase= ~Tx_phase;
}
#endif

//----------------------------------------------------------
#ifdef SWUART_RECEIVE_USED
/**
 * @brief Receive byte checking and store
 * @par Parameters:
 * u8 pointer to place where the received byte should be stored
 * @retval error code in low nibble (see codes definition in documentation 
 * - TRUE meens no error receive) or FALSE (no byte received)
 */
u8 uart_read(u8 *b) {
	u8 res;
	if(test_status(receive_buffer_full))	{
		*b= Rx_data;
#ifdef BIT9
		res= UART_sts & 0x2F;
#else
		res= UART_sts & 0x0F;
#endif
		UART_sts&= ~0x0F;
		return(res);
	}
	else
		return(FALSE);
}	
//----------------------------------------------------------
/**
 * @brief Timing of uart Rx - receive all bits proccess control
 * @par Parameters:
 * None
 * @retval None
 */
void uart_Rx_timing(void) {

}
#endif

/** 
  * @}
  * @}
  */
