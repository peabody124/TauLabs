/**
 ******************************************************************************
 * @addtogroup ESC esc
 * @brief The main ESC code
 *
 * @file       esc.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      INSGPS Test Program
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
#include "esc_fsm.h"
#include "fifo_buffer.h"
#include <pios_stm32.h>

#define CURRENT_LIMIT 1210

//TODO: Check the ADC buffer pointer and make sure it isn't dropping swaps
//TODO: Check the time commutation is being scheduled, make sure it's the future
//TODO: Slave two timers together so in phase
//TODO: Ideally lock ADC and delay timers together to both
//TODO: Look into using TIM1
//know the exact time of each sample and the PWM phase

//#define BACKBUFFER_ZCD
//#define BACKBUFFER_ADC
//#define BACKBUFFER_DIFF

/* Prototype of PIOS_Board_Init() function */
extern void PIOS_Board_Init(void);

#define DOWNSAMPLING 1

#if defined(BACKBUFFER_ADC) || defined(BACKBUFFER_ZCD) || defined(BACKBUFFER_DIFF)
uint16_t back_buf[2024];
uint16_t back_buf_point = 0;
#endif

#define LED_ERR LED1
#define LED_GO  LED2
#define LED_MSG LED3

int16_t zero_current = 0;

volatile uint8_t low_pin;
volatile uint8_t high_pin;
volatile uint8_t undriven_pin;
volatile bool pos;

const uint8_t dT = 1e6 / PIOS_ADC_RATE; // 6 uS per sample at 160k
float rate = 0;

static void test_esc();
static void panic(int diagnostic_code);
uint16_t pwm_duration ;
uint32_t counter = 0;

#define NUM_SETTLING_TIMES 20
uint32_t timer;
uint16_t timer_lower;
uint32_t step_period = 0x0080000;
uint32_t last_step = 0;
int16_t low_voltages[3];
int32_t avg_low_voltage;
struct esc_fsm_data * esc_data = 0;

/**
 * @brief ESC Main function
 */

uint32_t offs = 0;

int main()
{
	esc_data = 0;
	PIOS_Board_Init();


	PIOS_ADC_Config(DOWNSAMPLING);

	// TODO: Move this into a PIOS_DELAY function
	TIM_OCInitTypeDef tim_oc_init = {
		.TIM_OCMode = TIM_OCMode_PWM1,
		.TIM_OutputState = TIM_OutputState_Enable,
		.TIM_OutputNState = TIM_OutputNState_Disable,
		.TIM_Pulse = 0,
		.TIM_OCPolarity = TIM_OCPolarity_High,
		.TIM_OCNPolarity = TIM_OCPolarity_High,
		.TIM_OCIdleState = TIM_OCIdleState_Reset,
		.TIM_OCNIdleState = TIM_OCNIdleState_Reset,
	};
	TIM_OC1Init(TIM4, &tim_oc_init);
	TIM_ITConfig(TIM4, TIM_IT_CC1, ENABLE);  // Enabled by FSM

	NVIC_InitTypeDef NVIC_InitStructure;
	NVIC_InitStructure.NVIC_IRQChannel = TIM4_IRQn;
	NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = PIOS_IRQ_PRIO_HIGH;
	NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0;
	NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
	NVIC_Init(&NVIC_InitStructure);

	// This pull up all the ADC voltages so the BEMF when at -0.7V
	// is still positive
	PIOS_GPIO_Enable(0);
	PIOS_GPIO_Off(0);

	PIOS_LED_Off(LED_ERR);
	PIOS_LED_On(LED_GO);
	PIOS_LED_On(LED_MSG);

	test_esc();

	esc_data = esc_fsm_init();
	esc_data->speed_setpoint = 0;


	PIOS_ADC_StartDma();
	
	while(1) {
		counter++;

		esc_process_static_fsm_rxn();
	}
	return 0;
}

#define MAX_RUNNING_FILTER 512
uint32_t DEMAG_BLANKING = 20;

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
	74, // ... 1280
	68, // ... 1408
	62, // ... 1536
	57, // ... 1664
	53, // ... 1792
	50, // ... 1920
	47, // ... 2048
	44, // ... 2176
	41, // ... 2304
	39, // ... 2432
	37, // ... 2560
	35, // ... 2688
	34, // ... 2816
	32, // ... 2944
	31, // ... 3072
	30, // ... 3200
	29, // ... 3328
	33, // ... 3456
	32, // ... 3584
	31, // ... 3712
	30, // ... 3840
	29, // ... 3968
	28, // ... 4096
	27, // ... 4224
	26, // ... 4352
	25, // ... 4480
	24, // ... 4608
	23, // ... 4736
	22, // ... 4864
	21, // ... 4992
	21, // ... 5120
	20, // ... 5248
	19, // ... 5376
	18, // ... 5504
	18, // ... 5632
	18, // ... 5760
	18, // ... 5888
	17, // ... 6016
	17, // ... 6144
	17, // ... 6272
	17, // ... 6400
	17, // ... 6528
	16, // ... 6656
	16, // ... 6784
	16, // ... 6912
	16, // ... 7040
	8, // ... 7168
	8, // ... 7296
	8, // ... 7424
	8, // ... 7552
	7, // ... 7680
	7, // ... 7808
	7, // ... 7936
	5, // ... 8064
	5, // ... 8192
	5, // ... 8320
	5, // ... 8448
	5, // ... 8576
	5, // ... 8704
	5, // ... 8832
	5, // ... 8960
	4, // ... 9088
	4, // ... 9216
	4, // ... 9344
	4, // ... 9472
	4, // ... 9600
	4, // ... 9728
	4, // ... 9856
	4, // ... 9984
	4, // ... 10112
	4, // ... 10240
	4, // ... 10368
	4, // ... 10496
	4, // ... 10624
	4, // ... 10752
	4, // ... 10880
	4, // ... 11008
	4, // ... 11136
	4, // ... 11264
	4, // ... 11392
	4, // ... 11520
	4, // ... 11648
	3, // ... 11776
	3, // ... 11904
	2, // ... 12032
	2, // ... 12160
	2, // ... 12288
	2, // ... 12416
	2, // ... 12544
	2, // ... 12672
	2, // ... 12800
	2, // ... 12928
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

uint8_t filter_length_scalar = 100;

void DMA1_Channel1_IRQHandler(void)
{	
	static enum pios_esc_state prev_state = ESC_STATE_AB;
	static int16_t * raw_buf;
	enum pios_esc_state curr_state;

	if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->full_flag /*DMA1_IT_TC1*/)) {	// whole double buffer filled
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[pios_adc_devs[0].dma_half_buffer_size];
		DMA_ClearFlag(pios_adc_devs[0].cfg->full_flag);
//		PIOS_GPIO_On(1);		
	}
	else if (DMA_GetFlagStatus(pios_adc_devs[0].cfg->half_flag /*DMA1_IT_HT1*/)) {
		pios_adc_devs[0].valid_data_buffer = &pios_adc_devs[0].raw_data_buffer[0];
		DMA_ClearFlag(pios_adc_devs[0].cfg->half_flag);
//		PIOS_GPIO_Off(1); 
	}
	else {
		// This should not happen, probably due to transfer errors
		DMA_ClearFlag(pios_adc_devs[0].cfg->dma.irq.flags /*DMA1_FLAG_GL1*/);
		return;
	}

	bad_flips += (raw_buf == pios_adc_devs[0].valid_data_buffer);
	raw_buf = (int16_t *) pios_adc_devs[0].valid_data_buffer;

	curr_state = PIOS_ESC_GetState();

#ifdef BACKBUFFER_ADC
	// Debugging code - keep a buffer of old ADC values
	back_buf[back_buf_point++] = raw_buf[0];
	back_buf[back_buf_point++] = raw_buf[1];
	back_buf[back_buf_point++] = raw_buf[2];
	back_buf[back_buf_point++] = raw_buf[3];
	if(back_buf_point >= (NELEMENTS(back_buf)-3))
#endif

	if( PIOS_DELAY_DiffuS(esc_data->last_swap_time) < DEMAG_BLANKING )
		return;
	
	// Smooth the estimate of current a bit 	
	current_filter_sum += raw_buf[0];
	current_filter_sum -= current_filter[current_filter_pointer];
	current_filter[current_filter_pointer] = raw_buf[0];
	current_filter_pointer++;
	if(current_filter_pointer >= MAX_CURRENT_FILTER) 
		current_filter_pointer = 0;
#if MAX_CURRENT_FILTER == 16
	esc_data->current = (current_filter_sum >> 4) - zero_current;
#else
	esc_data->current = current_filter_sum / MAX_CURRENT_FILTER - zero_current;
#endif

	if(esc_data->current > CURRENT_LIMIT)
		esc_fsm_inject_event(ESC_EVENT_OVERCURRENT, 0);

	   
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
		
		if((esc_data->current_speed >> 7) > NELEMENTS(filter_length))
		   running_filter_length = filter_length[NELEMENTS(filter_length)-1];
		else
		   running_filter_length = filter_length[esc_data->current_speed >> 7];

		running_filter_length = running_filter_length * filter_length_scalar / 100;
		
		if(running_filter_length >= MAX_RUNNING_FILTER)
			running_filter_length = MAX_RUNNING_FILTER;
		
		if(esc_data->current_speed > 8000)
			running_filter_length = 10;

		switch(curr_state) {
			case ESC_STATE_AC:
				undriven_pin = 1;
				pos = true;
				break;
			case ESC_STATE_CA:
				undriven_pin = 1;
				pos = false;
				break;
			case ESC_STATE_AB:
				undriven_pin = 2;
				pos = false;
				break;
			case ESC_STATE_BA:
				undriven_pin = 2;
				pos = true;
				break;
			case ESC_STATE_BC:
				undriven_pin = 0;
				pos = false;
				break;
			case ESC_STATE_CB:
				undriven_pin = 0;
				pos = true;
				break;
			default:
				PIOS_ESC_Off();
		}
		init_time = PIOS_DELAY_DiffuS(timeval);
	}

	int32_t last_swap_delay = PIOS_DELAY_DiffuS(esc_data->last_swap_time);
	if(last_swap_delay >= 0 && last_swap_delay < DEMAG_BLANKING)
		return;

	calls_to_detect++;
	for(uint32_t i = 0; i < DOWNSAMPLING; i++) {
		int16_t undriven = raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + undriven_pin]; // - low_voltages[undriven_pin];
		int16_t ref = (raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + 0] + 
					   raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + 1] +
					   raw_buf[PIOS_ADC_NUM_CHANNELS * i + 1 + 2] - avg_low_voltage*0) / 3;
		int32_t diff;
		
		if(pos)
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
#ifdef BACKBUFFER_ZCD
			back_buf[back_buf_point++] = below_time;
			back_buf[back_buf_point++] = esc_data->speed_setpoint;
			if(back_buf_point > (sizeof(back_buf) / sizeof(back_buf[0])))
				back_buf_point = 0;
#endif /* BACKBUFFER_ZCD */
			break;
		}
#ifdef BACKBUFFER_DIFF
		// Debugging code - keep a buffer of old ADC values
		back_buf[back_buf_point++] = diff;
		if(back_buf_point >= NELEMENTS(back_buf))
			back_buf_point = 0;
#endif
	}	
}

/* INS functions */
void panic(int diagnostic_code)
{
	PIOS_ESC_Off();
	// Polarity backwards
	PIOS_LED_On(LED_ERR);
	while(1) {
		for(int i=0; i<diagnostic_code; i++)
		{
			PIOS_LED_Toggle(LED_ERR);
			PIOS_DELAY_WaitmS(250);
			PIOS_LED_Toggle(LED_ERR);
			PIOS_DELAY_WaitmS(250);
		}
		PIOS_DELAY_WaitmS(1000);
	}
}

//TODO: Abstract out constants.  Need to know battery voltage too
int32_t voltages[6][3];
void test_esc() {


	PIOS_ESC_Off();
	PIOS_DELAY_WaitmS(150);
	zero_current = PIOS_ADC_PinGet(0);

	PIOS_ESC_Arm();

	PIOS_ESC_TestGate(ESC_A_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[0] = PIOS_ADC_PinGet(1);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[1] = PIOS_ADC_PinGet(2);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	low_voltages[2] = PIOS_ADC_PinGet(3);
	avg_low_voltage = low_voltages[0] + low_voltages[1] + low_voltages[2];

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_A_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[1][0] = PIOS_ADC_PinGet(1);
	voltages[1][1] = PIOS_ADC_PinGet(2);
	voltages[1][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_A_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[0][0] = PIOS_ADC_PinGet(1);
	voltages[0][1] = PIOS_ADC_PinGet(2);
	voltages[0][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[3][0] = PIOS_ADC_PinGet(1);
	voltages[3][1] = PIOS_ADC_PinGet(2);
	voltages[3][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_B_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[2][0] = PIOS_ADC_PinGet(1);
	voltages[2][1] = PIOS_ADC_PinGet(2);
	voltages[2][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_LOW);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[5][0] = PIOS_ADC_PinGet(1);
	voltages[5][1] = PIOS_ADC_PinGet(2);
	voltages[5][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE / 2);
	PIOS_ESC_TestGate(ESC_C_HIGH);
	PIOS_DELAY_WaituS(250);
	PIOS_ESC_SetDutyCycle(PIOS_ESC_MAX_DUTYCYCLE);
	PIOS_DELAY_WaituS(300);
	voltages[4][0] = PIOS_ADC_PinGet(1);
	voltages[4][1] = PIOS_ADC_PinGet(2);
	voltages[4][2] = PIOS_ADC_PinGet(3);

	PIOS_ESC_Off();
	// If the particular phase isn't moving fet is dead
	if(voltages[0][0] < 1000)
		panic(1);
	if(voltages[1][0] > 900)
		panic(2);
	if(voltages[2][1] < 1000)
		panic(2);
	if(voltages[3][1] > 900)
		panic(3);
	if(voltages[4][2] < 1000)
		panic(4);
	if(voltages[5][2] > 900)
		panic(5);

	// TODO: If other channels don't follow then motor lead bad
}

uint32_t bad_inputs;
void PIOS_TIM_4_irq_override();
extern void PIOS_DELAY_timeout();
void TIM4_IRQHandler(void) __attribute__ ((alias ("PIOS_TIM_4_irq_handler")));
static bool rising = false;

static uint16_t rise_value;
static uint16_t fall_value;
static uint16_t capture_value;

static void PIOS_TIM_4_irq_handler (void)
{
	static uint32_t last_input_update;
	
	if(TIM_GetITStatus(TIM4,TIM_IT_CC1)) {
		PIOS_DELAY_timeout();
		TIM_ClearITPendingBit(TIM4,TIM_IT_CC1);
	}
	
	if (TIM_GetITStatus(TIM4, TIM_IT_CC3)) {
		
		TIM_ClearITPendingBit(TIM4,TIM_IT_CC3);
		
		TIM_ICInitTypeDef TIM_ICInitStructure = {
			.TIM_ICPolarity = TIM_ICPolarity_Rising,
			.TIM_ICSelection = TIM_ICSelection_DirectTI,
			.TIM_ICPrescaler = TIM_ICPSC_DIV1,
			.TIM_ICFilter = 0x0,
		};
		
		uint16_t tmp;
		if(rising) {
			rising = false;
			rise_value = TIM_GetCapture3(TIM4);
			
			tmp = (rise_value & 0xFFFC) | 0x0001;
			PIOS_COM_SendBuffer(PIOS_COM_DEBUG, (uint8_t *) &tmp, 2);

			/* Switch polarity of input capture */
			TIM_ICInitStructure.TIM_ICPolarity = TIM_ICPolarity_Falling;
			TIM_ICInitStructure.TIM_Channel = TIM_Channel_3;
			TIM_ICInit(TIM4, &TIM_ICInitStructure);
		} else {
			rising = true;
			fall_value = TIM_GetCapture3(TIM4);

			tmp = (fall_value & 0xFFFC) | 0x0002;
			PIOS_COM_SendBuffer(PIOS_COM_DEBUG, (uint8_t *) &tmp, 2);

			/* Switch polarity of input capture */
			TIM_ICInitStructure.TIM_ICPolarity = TIM_ICPolarity_Rising;
			TIM_ICInitStructure.TIM_Channel = TIM_Channel_3;
			TIM_ICInit(TIM4, &TIM_ICInitStructure);
			
			if (fall_value > rise_value) {
				capture_value = fall_value - rise_value;
			} else {
				capture_value = TIM4->ARR + fall_value - rise_value;
			}
			
			capture_value &= 0xFFFC;
			PIOS_COM_SendBuffer(PIOS_COM_DEBUG, (uint8_t *) &capture_value, 2);

			last_input_update = PIOS_DELAY_GetRaw();
		}
		
		
		if(capture_value < 900 || capture_value > 2200)
			capture_value = 0;
		else {
			esc_data->speed_setpoint = (capture_value < 1050) ? 0 : 400 + (capture_value - 1050) * 6;
		}
	} 
	
	if (TIM_GetITStatus(TIM4, TIM_IT_Update)) {
		if (PIOS_DELAY_DiffuS(last_input_update) > 100000)
			esc_data->speed_setpoint = 0;
		TIM_ClearITPendingBit(TIM4,TIM_IT_Update);
	}
}

/*
 Notes:
 1. For start up, definitely want to use complimentary PWM to ground the lower side, making zero crossing truly "zero"
 2. May want to use the "middle" sensor to actually pull it up, so that zero is above zero (in ADC range).  Should still
    see BEMF at -0.7 (capped by transistor range) relative to that point (divided down by whatever)
 3. Possibly use an inadequate voltage divider plus use the TVS cap to keep the part of the signal near zero clean
 */


