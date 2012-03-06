/**
 ******************************************************************************
 * @addtogroup ESC esc
 * @brief The ESC ADC processing
 *
 * @file       esc_adc.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2012.
 * @brief      Code to process the ADC data to detect BEMF crossing
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

/* OpenPilot Includes */
#include "pios.h"
#include "esc.h"

//! External global variables about the ESC information
extern struct esc_fsm_data * esc_data;
extern struct esc_control esc_control;

#define MAX_RUNNING_FILTER 512
uint32_t DEMAG_BLANKING = 70;

#define MAX_CURRENT_FILTER 16
int32_t current_filter[MAX_CURRENT_FILTER];
int32_t current_filter_sum = 0;
int32_t current_filter_pointer = 0;

static int16_t diff_filter[MAX_RUNNING_FILTER];
static int32_t diff_filter_pointer = 0;
int32_t running_filter_length = 512;
static int32_t running_filter_sum = 0;

uint32_t calls_to_detect = 0;
uint32_t calls_to_last_detect = 0;

uint32_t samples_averaged;
int32_t threshold = -10;
int32_t divisor = 40;
#include "pios_adc_priv.h"
uint32_t detected;
uint32_t bad_flips;

uint32_t init_time;

// Index to look up filter lengths by RPM (incrementing by 128 rpm per entry
static const int32_t filter_length[] = {
	744, // ... 128
	372, // ... 256
	248, // ... 384
	186, // ... 512
	149, // ... 640
	124, // ... 768
	106, // ... 896
	93, // ... 1024
	83, // ... 1152
	64, // ... 1280
	58, // ... 1408
	52, // ... 1536
	47, // ... 1664
	43, // ... 1792
	40, // ... 1920
	37, // ... 2048
	34, // ... 2176
	31, // ... 2304
	29, // ... 2432
	27, // ... 2560
	25, // ... 2688
	24, // ... 2816
	22, // ... 2944
	21, // ... 3072
	20, // ... 3200
	19, // ... 3328
	18, // ... 3456
	18, // ... 3584
	17, // ... 3712
	17, // ... 3840
	17, // ... 3968
	16, // ... 4096
	16, // ... 4224
	16, // ... 4352
	15, // ... 4480
	15, // ... 4608
	15, // ... 4736
	14, // ... 4864
	14, // ... 4992
	14, // ... 5120
	14, // ... 5248
	14, // ... 5376
	13, // ... 5504
	13, // ... 5632
	13, // ... 5760
	13, // ... 5888
	13, // ... 6016
	12, // ... 6144
	12, // ... 6272
	12, // ... 6400
	12, // ... 6528
	12, // ... 6656
	12, // ... 6784
	12, // ... 6912
	11, // ... 7040
	11, // ... 7168
	11, // ... 7296
	11, // ... 7424
	11, // ... 7552
	10, // ... 7680
	10, // ... 7808
	10, // ... 7936
	9, // ... 8064
	9, // ... 8192
	9, // ... 8320
	9, // ... 8448
	9, // ... 8576
	9, // ... 8704
	9, // ... 8832
	8, // ... 8960
	8, // ... 9088
	8, // ... 9216
	8, // ... 9344
	7, // ... 9472
	7, // ... 9600
	7, // ... 9728
	7, // ... 9856
	6, // ... 9984
	6, // ... 10112
	6, // ... 10240
	6, // ... 10368
	6, // ... 10496
	6, // ... 10624
	6, // ... 10752
	5, // ... 10880
	5, // ... 11008
	5, // ... 11136
	5, // ... 11264
	5, // ... 11392
	5, // ... 11520
	5, // ... 11648
	4, // ... 11776
	4, // ... 11904
	4, // ... 12032
	4, // ... 12160
	4, // ... 12288
	4, // ... 12416
	3, // ... 12544
	3, // ... 12672
	3, // ... 12800
	3, // ... 12928
	2, // ... 13056
	2, // ... 13184
	2, // ... 13312
	2, // ... 13440
	2, // ... 13568
	2, // ... 13696
	2, // ... 13824
	2, // ... 13952
	2, // ... 14080
	2, // ... 14208
	2, // ... 14336
	2, // ... 14464
	2, // ... 14592
	2, // ... 14720
	2, // ... 14848
	2, // ... 14976
	2, // ... 15104
	2, // ... 15232
	2, // ... 15360
	2, // ... 15488
	2, // ... 15616
	2, // ... 15744
	2, // ... 15872
	2, // ... 16000
	2, // ... 16128
	2, // ... 16256
	2, // ... 16384
	2, // ... 16512
	2, // ... 16640
	2, // ... 16768
	2, // ... 16896
	2, // ... 17024
	2, // ... 17152
	2, // ... 17280
	2, // ... 17408
	2, // ... 17536
	2, // ... 17664
	2, // ... 17792
	2, // ... 17920
	1,
};	

uint8_t filter_length_scalar = 25;

const uint8_t undriven_pin[6] = {
	[ESC_STATE_AC] = 1,
	[ESC_STATE_CA] = 1,
	[ESC_STATE_AB] = 2,
	[ESC_STATE_BA] = 2,
	[ESC_STATE_BC] = 0,
	[ESC_STATE_CB] = 0
};

// This depends on the driving mode, assuming PWM high
const uint8_t high_pin[6] = {
	[ESC_STATE_AC] = 2,
	[ESC_STATE_CA] = 0,
	[ESC_STATE_AB] = 1,
	[ESC_STATE_BA] = 0,
	[ESC_STATE_BC] = 2,
	[ESC_STATE_CB] = 1
};

const bool positive_forward[6] = {
	[ESC_STATE_AC] = true,
	[ESC_STATE_CA] = false,
	[ESC_STATE_AB] = false,
	[ESC_STATE_BA] = true,
	[ESC_STATE_BC] = false,
	[ESC_STATE_CB] = true
};

const bool positive_backwards[6] = {
	[ESC_STATE_AC] = false,
	[ESC_STATE_CA] = true,
	[ESC_STATE_AB] = true,
	[ESC_STATE_BA] = false,
	[ESC_STATE_BC] = true,
	[ESC_STATE_CB] = false
};

uint32_t adc_count = 0;
void DMA1_Channel1_IRQHandler(void)
{	
	adc_count++;	
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	//	static uint8_t overcurrent_count = 0;
	static int16_t * raw_buf;
	enum pios_esc_state curr_state;
	
	if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->full_flag /*DMA1_IT_TC1*/)) {	// whole double buffer filled
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[pios_adc_devs[0].dma_half_buffer_size];
		DMA_ClearFlag(pios_adc_devs[0].cfg->full_flag);
	}
	else if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->half_flag /*DMA1_IT_HT1*/)) {
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[0];
		DMA_ClearFlag(pios_adc_devs[0].cfg->half_flag);
	}
	else {
		// This should not happen, probably due to transfer errors
		DMA_ClearFlag(pios_adc_devs[0].cfg->dma.irq.flags /*DMA1_FLAG_GL1*/);
		return;
	}
	
	if(!esc_data)
		return;
	
	bad_flips += (raw_buf == pios_adc_devs[0].valid_data_buffer);
	raw_buf = (int16_t *) pios_adc_devs[0].valid_data_buffer;
	
	curr_state = PIOS_ESC_GetState();
	
	// If requested log data until buffer full
	if(esc_control.backbuffer_logging_status == ESC_LOGGING_CAPTURE) {
		// When T3 flag is off then it is sampled when driven
		uint16_t flags = 0;
		flags |= (TIM2->CR1 & TIM_CR1_DIR) ? 0x8000 : 0;
		flags |= (TIM3->CR1 & TIM_CR1_DIR) ? 0x4000 : 0;
		int16_t samples[4] = {raw_buf[0], raw_buf[1] | flags, raw_buf[2] | flags, raw_buf[3] | flags};
		if (esc_logger_put((uint8_t *) &samples[0], 8) != 0)
			esc_control.backbuffer_logging_status = ESC_LOGGING_FULL;
	}
	
	if( PIOS_DELAY_DiffuS(esc_data->last_swap_time) < DEMAG_BLANKING )
		return;
	
	// Smooth the estimate of current a bit 	
	int32_t current_measurement = raw_buf[0] * 8;// Convert to almost mA
	esc_data->current_ma = (127 * esc_data->current_ma + current_measurement) / 128; // IIR filter
	
	/*	if (esc_data->current > CURRENT_LIMIT && overcurrent_count > 4)
	 esc_fsm_inject_event(ESC_EVENT_OVERCURRENT, 0);
	 else if (esc_data->current > CURRENT_LIMIT)
	 overcurrent_count ++;
	 else overcurrent_count = 0; */
	
	
	// If detected this commutation don't bother here
	if(curr_state == prev_state && esc_data->detected)
		return;
	else if(curr_state != prev_state) {
		uint32_t timeval = PIOS_DELAY_GetRaw();
		
		prev_state = curr_state;
		calls_to_detect = 0;
		running_filter_sum = 0;
		diff_filter_pointer = 0;
		samples_averaged = 0;
		
		// Currently can get into a situation at high RPMs that causes filter to get
		// too short and then trigger lots of false ZCDs when it thinks its running at
		// 15k rpm
		uint32_t current_speed = esc_data->current_speed;
		if (current_speed > 8500)
			current_speed = 8500;
		
		if((current_speed >> 7) > NELEMENTS(filter_length))
			running_filter_length = filter_length[NELEMENTS(filter_length)-1];
		else
			running_filter_length = filter_length[esc_data->current_speed >> 7];
		
		running_filter_length = running_filter_length * filter_length_scalar / 100;
		
		if(running_filter_length >= MAX_RUNNING_FILTER)
			running_filter_length = MAX_RUNNING_FILTER;
		
		init_time = PIOS_DELAY_DiffuS(timeval);
	}
	
	int32_t last_swap_delay = PIOS_DELAY_DiffuS(esc_data->last_swap_time);
	if(last_swap_delay >= 0 && last_swap_delay < DEMAG_BLANKING)
		return;
	
	calls_to_detect++;
	
	int16_t undriven = raw_buf[1 + undriven_pin[curr_state]];
	int16_t high = raw_buf[1 + high_pin[curr_state]];
	int16_t ref = (raw_buf[1 + 0] + 
				   raw_buf[1 + 1] +
				   raw_buf[1 + 2]) / 3;
	int32_t diff;
	
	if(esc_data->duty_cycle > (PIOS_ESC_MAX_DUTYCYCLE * 0.05)) {
		// 127 / 27 is the 12.7 / 2.7 resistor divider
		// 3300 / 4096 is the mv / LSB
		// (high * PIOS_ESC_MAX_DUTYCYCLE * 127 * 3300) / (esc_data->duty_cycle * 27 * 4096);
		// Shifted equation around to not run out of precision
		uint32_t battery_mv = (high * 127 * 3300 / (esc_data->duty_cycle * 27)) / (4096 / PIOS_ESC_MAX_DUTYCYCLE);
		esc_data->battery_mv = (16383 * esc_data->battery_mv + battery_mv) / 16384;
	}
	
	const bool *pos = (PIOS_ESC_GetDirection() == ESC_FORWARD) ? positive_forward : positive_backwards;
	if(pos[curr_state])
		diff = undriven - ref;
	else
		diff = ref - undriven;
	
	// Update running sum and history
	running_filter_sum += diff;
	// To avoid having to wipe the filter each commutation
	if(samples_averaged >= running_filter_length)
		running_filter_sum -= diff_filter[diff_filter_pointer];
	diff_filter[diff_filter_pointer] = diff;
	diff_filter_pointer++;
	if(diff_filter_pointer >= running_filter_length)
		diff_filter_pointer = 0;
	samples_averaged++;
	
	if(running_filter_sum > 0 && samples_averaged >= running_filter_length) {
		esc_data->detected = true;
		esc_fsm_inject_event(ESC_EVENT_ZCD, TIM4->CNT);
		calls_to_last_detect = calls_to_detect;
	}	
}

