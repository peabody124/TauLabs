/*
 *  esc_fsm.c
 *  OpenPilotOSX
 *
 *  Created by James Cotton on 6/12/11.
 *  Copyright 2011 OpenPilot. All rights reserved.
 *
 */

#include "esc_fsm.h"
#include "pios.h"
#include "stdbool.h"
#include "pios_esc.h"

// TOOD: Important safety features
// 1. Disable after a timeout from last control command
// 2. Reliable shutdown on high current
// 3. Monitor high voltage level to track battery (bonus for determining # of cells)
// 4. BEEPS!

// 5. In case of soft current limit, go into current limiting mode
// 6. Freehweeling mode for when a ZCD is missed

#define RPM_TO_US(x) (1e6 * 60 / (x * COMMUTATIONS_PER_ROT) )
#define US_TO_RPM(x) RPM_TO_US(x)
#define COMMUTATIONS_PER_ROT (7*6)
#define ESC_CONFIG_MAGIC 0x763fedc

struct esc_config config = {
	.kp = 0.0025,
	.ki = 0.001,
	.kff = 1.3e-4,
	.kff2 = -0.05,
	.ilim = 0.5,
	.max_dc_change = 20,
	.min_dc = 0.01,
	.max_dc = 0.99,
	.initial_startup_speed = 400,
	.final_startup_speed = 1200,
	.commutation_phase = 0.5,
	.soft_current_limit = 250,
	.hard_current_limit = 850,
	.magic = ESC_CONFIG_MAGIC,
};

static void go_esc_nothing(uint16_t time) {};
// Fault and stopping transitions
static void go_esc_fault(uint16_t);
static void go_esc_stopping(uint16_t);
static void go_esc_stopped(uint16_t);

// Startup states
static void go_esc_startup_enable(uint16_t);
static void go_esc_startup_grab(uint16_t);
static void go_esc_startup_wait(uint16_t);
static void go_esc_startup_zcd(uint16_t);
static void go_esc_startup_nozcd(uint16_t);

// Closed-loop states
static void go_esc_cl_start(uint16_t time);
static void go_esc_cl_commutated(uint16_t time);
static void go_esc_cl_zcd(uint16_t time);
static void go_esc_cl_nozcd(uint16_t time);


// Misc utility functions
static void esc_fsm_process_auto(uint16_t);
static void esc_fsm_schedule_event(enum esc_event event, uint16_t time);

static void commutate();
static void zcd(uint16_t time);

const static struct esc_transition esc_transition[ESC_FSM_NUM_STATES] = {
	[ESC_STATE_FSM_FAULT] = {
		.entry_fn = go_esc_fault,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_INIT] = {
		.entry_fn = go_esc_nothing,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_WAIT_FOR_ARM,
		},
	},
	[ESC_STATE_WAIT_FOR_ARM] = {
		.entry_fn = go_esc_nothing,
		.next_state = {
			[ESC_EVENT_ARM] = ESC_STATE_IDLE,
		},
	},
	[ESC_STATE_IDLE] = {
		.entry_fn = go_esc_nothing,
		.next_state = {
			[ESC_EVENT_START] = ESC_STATE_STARTUP_ENABLE,
		},
	},
	[ESC_STATE_STOPPING] = {
		.entry_fn = go_esc_stopping,
		.next_state = {
			// TODO: Find interesting criterion to consider stopped on
			[ESC_EVENT_AUTO] = ESC_STATE_STOPPED,
		},
	},
	[ESC_STATE_STOPPED] = {
		.entry_fn = go_esc_stopped,
		.next_state = {
			// For now just loop in stopped state while we auto running
			[ESC_EVENT_AUTO] = ESC_STATE_STOPPED /*ESC_STATE_IDLE*/,
		},
	},

	/* Startup states */
	[ESC_STATE_STARTUP_ENABLE] = {
		.entry_fn = go_esc_startup_enable,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_STARTUP_GRAB,
		},
	},
	[ESC_STATE_STARTUP_GRAB] = {
		.entry_fn = go_esc_startup_grab,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_WAIT,
		},
	},
	[ESC_STATE_STARTUP_WAIT] = {
		// Schedule commutation
		.entry_fn = go_esc_startup_wait,
		.next_state = {
			// Could not get BEMF in time
			[ESC_EVENT_TIMEOUT] = ESC_STATE_STOPPING,
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_NOZCD_COMMUTATED,
		},
	},
	[ESC_STATE_STARTUP_ZCD_DETECTED] = {
		.entry_fn = go_esc_startup_zcd,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_WAIT,
			[ESC_EVENT_CLOSED] = ESC_STATE_CL_START,
		},
	},
	[ESC_STATE_STARTUP_NOZCD_COMMUTATED] = {
		.entry_fn = go_esc_startup_nozcd,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_NOZCD_COMMUTATED,
		},
	},

	/* Closed loop states */
	[ESC_STATE_CL_START] = {
		.entry_fn = go_esc_cl_start,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},
	[ESC_STATE_CL_COMMUTATED] = {
		.entry_fn = go_esc_cl_commutated,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_CL_ZCD,
			[ESC_EVENT_TIMEOUT] = ESC_STATE_CL_NOZCD,
		},
	},
	[ESC_STATE_CL_NOZCD] = {
		.entry_fn = go_esc_cl_nozcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},
	[ESC_STATE_CL_ZCD] = {
		.entry_fn = go_esc_cl_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},

};

static struct esc_fsm_data esc_data;

void esc_process_static_fsm_rxn() {

	if(esc_data.current > config.hard_current_limit)
		esc_fsm_inject_event(ESC_EVENT_FAULT, 0);

	switch(esc_data.state) {
		case ESC_STATE_FSM_FAULT:
			// Shouldn't happen
			break;
		case ESC_STATE_INIT:
			// No static rxn.  Goes straight to wait_for_arm.
			break;
		case ESC_STATE_WAIT_FOR_ARM:
			// TODO: Monitor for low input.
			// Now now make it immediately arm
			esc_fsm_inject_event(ESC_EVENT_ARM,0);
			break;
		case ESC_STATE_IDLE:
			// TODO: Monitor for non-low input
			// For now just start
			esc_fsm_inject_event(ESC_EVENT_START,0);
			break;
		case ESC_STATE_STOPPING:
			// Monitor that the motor has stopped and current < threshold
			PIOS_ESC_Off();
			break;
		case ESC_STATE_STOPPED:
			// Not used
			break;

		/**
		 * Startup state transitions
		 */
		case ESC_STATE_STARTUP_GRAB:
		case ESC_STATE_STARTUP_WAIT:
		case ESC_STATE_STARTUP_ZCD_DETECTED:
		case ESC_STATE_STARTUP_NOZCD_COMMUTATED:
		{
			// Timing adjusted in entry function
			// This should perform run a current control loop
			float current_error = (100 - esc_data.current);
			current_error *= 0.0001;
			if(current_error > 0.0002)
				current_error = 0.0002;
			if(current_error < -0.0002)
				current_error = -0.0002;
			esc_data.duty_cycle += current_error;

			if(esc_data.duty_cycle < 0.01)
				esc_data.duty_cycle = 0.01;
			if(esc_data.duty_cycle > 0.2)
				esc_data.duty_cycle = 0.2;
			PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
			break;
		}

		/**
		 * Closed loop transitions
		 */
		case ESC_STATE_CL_START:
		case ESC_STATE_CL_COMMUTATED:
		case ESC_STATE_CL_NOZCD:
		case ESC_STATE_CL_ZCD:
		{
			static uint16_t last_timer;
			uint16_t cur_timer = PIOS_DELAY_GetuS();

			if(esc_data.current > config.soft_current_limit) {

				float dT = (uint16_t) (cur_timer - last_timer);
				dT *= 1e-6; // convert to seconds
				float max_dc_change = config.max_dc_change * dT;

				esc_data.duty_cycle -= max_dc_change;
				PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
			}
			last_timer = cur_timer;
		}
			break;

		// By default do nothing
		default:
			break;
	}
}

/**
 * Invalid state transition occurred.
 */
static void go_esc_fault(uint16_t time)
{
	// Stop ADC quickly for grabbing debug buffer
	PIOS_ADC_StopDma();
	PIOS_ESC_Off();
	esc_data.faults ++;
}

/**
 * Start shutting down the motors
 */
static void go_esc_stopping(uint16_t time)
{
	PIOS_ESC_Off();
}

/**
 * Motors successfully stopped
 */
static void go_esc_stopped(uint16_t time)
{
	TIM_ITConfig(TIM4, TIM_IT_CC1, DISABLE);
	PIOS_ADC_StopDma();
}

/**
 * Enable ADC, timers, drivers, etc
 */
static void go_esc_startup_enable(uint16_t time)
{
	PIOS_ADC_StartDma();
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);

	// TODO: This should not have hardware calls
	TIM_ITConfig(TIM4, TIM_IT_CC1, ENABLE);

	PIOS_ESC_Arm();
}

/**
 * Grab the stator for startup sequencing
 */
static void go_esc_startup_grab(uint16_t time)
{
	// TODO: Set up a timeout for whole startup system

	PIOS_ESC_SetState(0);
	esc_data.current_speed = 200;
	esc_data.duty_cycle = 0.10;
	PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, PIOS_DELAY_GetuS() + 30000);  // Grab stator for 30 ms
}

/**
 * Called each time we enter the startup wait state or when a commutation occurs after ZCD
 * schedule a commutation.
 */
static void go_esc_startup_wait(uint16_t time)
{
	commutate();
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, PIOS_DELAY_GetuS() + RPM_TO_US(esc_data.current_speed));
}

/**
 * ZCD detected. Hang out here till commutate.
 */
static void go_esc_startup_zcd(uint16_t time)
{
	esc_data.consecutive_detected++;
	esc_data.consecutive_missed = 0;

	zcd(time);

	// Since we aren't getting ZCD keep accelerating
	if(esc_data.current_speed < config.final_startup_speed)
		esc_data.current_speed+=0;

	// Schedule next commutation
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, time + RPM_TO_US(esc_data.current_speed) / 2);

	if(esc_data.consecutive_detected > 20) {
		esc_fsm_inject_event(ESC_EVENT_CLOSED, time);
	}
}

/**
 * Commutation occured during startup without a ZCD
 */
static void go_esc_startup_nozcd(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed++;

	if(esc_data.consecutive_missed > 10) {
		esc_fsm_inject_event(ESC_EVENT_FAULT, time);
	} else {
		commutate();

		// Since we aren't getting ZCD keep accelerating
		if(esc_data.current_speed < config.final_startup_speed)
			esc_data.current_speed+=10;

		// Schedule next commutation
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, PIOS_DELAY_GetuS() + RPM_TO_US(esc_data.current_speed));
	}
}

/**
 * Commutation in closed loop
 */
static void go_esc_cl_start(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed = 0;
}

/**
 * Commutation in closed loop
 */
static void go_esc_cl_commutated(uint16_t time)
{
	commutate();
	esc_fsm_schedule_event(ESC_EVENT_TIMEOUT, time + RPM_TO_US(esc_data.current_speed) * 2);
	esc_data.Kv = esc_data.current_speed / esc_data.current;

//	if(esc_data.Kv < 15)
//		esc_fsm_inject_event(ESC_EVENT_FAULT, 0);
}

/**
 * When a zcd is detected
 */
static void go_esc_cl_zcd(uint16_t time)
{
	esc_data.consecutive_detected++;
	esc_data.consecutive_missed = 0;

	esc_data.dT = (uint16_t) (time - esc_data.last_zcd_time);
	esc_data.dT *= 1e-6; // convert to seconds

	float max_dc_change = config.max_dc_change * esc_data.dT;

	if(esc_data.current > config.soft_current_limit) {
		esc_data.duty_cycle -= max_dc_change;
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	} else {
		uint32_t interval = 0;
		for (uint8_t i = 0; i < NUM_STORED_SWAP_INTERVALS; i++)
			interval += esc_data.swap_intervals[i];
		interval /= NUM_STORED_SWAP_INTERVALS;
		esc_data.current_speed = US_TO_RPM(interval);

		int16_t error = esc_data.speed_setpoint - US_TO_RPM(interval);

		esc_data.error_accum += config.ki * error * esc_data.dT;
		if(esc_data.error_accum > config.ilim)
			esc_data.error_accum = config.ilim;
		if(esc_data.error_accum < -config.ilim)
			esc_data.error_accum = -config.ilim;

		float new_dc = esc_data.speed_setpoint * config.kff - config.kff2 + config.kp * error + esc_data.error_accum;

		if((new_dc - esc_data.duty_cycle) > max_dc_change)
			esc_data.duty_cycle += max_dc_change;
		else if((new_dc - esc_data.duty_cycle) < -max_dc_change)
			esc_data.duty_cycle -= max_dc_change;
		else
			esc_data.duty_cycle = new_dc;

		if(esc_data.duty_cycle < config.min_dc)
			esc_data.duty_cycle = config.min_dc;
		if(esc_data.duty_cycle > config.max_dc)
			esc_data.duty_cycle = config.max_dc;
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);

		zcd(time);

	}

	uint16_t zcd_delay = RPM_TO_US(esc_data.current_speed) * (1 - config.commutation_phase);
	int32_t future = time + zcd_delay - PIOS_DELAY_GetuS();
	if(future < 10)
		esc_fsm_inject_event(ESC_EVENT_FAULT, 0);
	else
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, time + (zcd_delay));
}

/**
 * For now fault out if we miss a zcd
 */
static void go_esc_cl_nozcd(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed++;
	if(esc_data.consecutive_missed > 10)
		esc_fsm_inject_event(ESC_EVENT_FAULT, 0);
	else
		esc_fsm_inject_event(ESC_EVENT_COMMUTATED, 0);
	esc_data.total_missed++;

}

/**
 * Initialize the FSM
 */
struct esc_fsm_data * esc_fsm_init()
{
	PIOS_ESC_Off();
	esc_data.state = ESC_STATE_INIT;
	esc_data.duty_cycle = 0;
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed = 0;
	esc_data.current_speed = 0;
	for (uint8_t i = 0; i < NUM_STORED_SWAP_INTERVALS; i++) {
		esc_data.swap_intervals[i] = 0;
		esc_data.zcd_intervals[i] = 0;
	}
	esc_data.swap_intervals_pointer = 0;
	esc_data.zcd_intervals_pointer = 0;
	esc_data.last_swap_time = 0;
	esc_data.last_zcd_time = 0;
	esc_data.scheduled_event_armed = false;
	esc_data.faults = 0;
	esc_data.total_missed = 0;

	esc_fsm_process_auto(PIOS_DELAY_GetuS());

	return &esc_data;
}


/**
 * Respond to an event with state system
 */
void esc_fsm_inject_event(enum esc_event event, uint16_t time)
{
	PIOS_IRQ_Disable();

	/*
	 * Move to the next state
	 *
	 * This is done prior to calling the new state's entry function to
	 * guarantee that the entry function never depends on the previous
	 * state.  This way, it cannot ever know what the previous state was.
	 */
	esc_data.state = esc_transition[esc_data.state].next_state[event];

	/* Call the entry function (if any) for the next state. */
	if (esc_transition[esc_data.state].entry_fn) {
		esc_transition[esc_data.state].entry_fn(time);
	}

	/* Process any AUTO transitions in the FSM */
	esc_fsm_process_auto(time);

	PIOS_IRQ_Enable();
}

/**
 * Process any automatic transitions
 */
static void esc_fsm_process_auto(uint16_t time)
{
	PIOS_IRQ_Disable();

	while (esc_transition[esc_data.state].next_state[ESC_EVENT_AUTO]) {
		esc_data.state = esc_transition[esc_data.state].next_state[ESC_EVENT_AUTO];

		/* Call the entry function (if any) for the next state. */
		if (esc_transition[esc_data.state].entry_fn) {
			esc_transition[esc_data.state].entry_fn(time);
		}
	}

	PIOS_IRQ_Enable();
}


//////////////////////////////////////////////////////////////////
//   More hardware dependent than I'd like
//   General functions not really for FSM explicitly
//////////////////////////////////////////////////////////////////
// TODO: Decouple these two functions from TIM4

/**
 * Schedule an interrupt from the timer to inject an event
 */
static void esc_fsm_schedule_event(enum esc_event event, uint16_t time)
{
	PIOS_IRQ_Disable();
	esc_data.scheduled_event = event;
	esc_data.scheduled_event_armed = true;
	TIM_SetCompare1(TIM4, time);
	PIOS_IRQ_Enable();
}

/**
 * The event interrupt
 */
void PIOS_DELAY_timeout() {
	TIM_ClearITPendingBit(TIM4,TIM_IT_CC1);
	TIM_ClearFlag(TIM4,TIM_IT_CC1);
	if(esc_data.scheduled_event_armed) {
		esc_data.scheduled_event_armed = false;
		esc_fsm_inject_event(esc_data.scheduled_event, PIOS_DELAY_GetuS());
	}
}

/**
 * Function to commutate the ESC and store any necessary information
 */
static void commutate()
{
	//PIOS_COM_SendFormattedStringNonBlocking(PIOS_COM_DEBUG, "%u %u\n", (next - swap_time), (uint32_t) (dc * 10000));
	uint16_t this_swap_time = PIOS_DELAY_GetuS();
	esc_data.last_swap_interval = (uint16_t) (this_swap_time - esc_data.last_swap_time);
	esc_data.last_swap_time = this_swap_time;

	esc_data.swap_intervals[esc_data.swap_intervals_pointer++] = esc_data.last_swap_interval;
	if(esc_data.swap_intervals_pointer > NUM_STORED_SWAP_INTERVALS)
		esc_data.swap_intervals_pointer = 0;

	esc_data.detected = false;

	PIOS_ESC_NextState();
}

/**
 * Function to handle the book-keeping when a zcd detected
 */
static void zcd(uint16_t time)
{
	uint16_t swap_diff = (uint16_t) (time - esc_data.last_swap_time);
	esc_data.last_zcd_time = time;

	esc_data.zcd_intervals[esc_data.zcd_intervals_pointer++] = swap_diff;
	if(esc_data.zcd_intervals_pointer > NUM_STORED_SWAP_INTERVALS)
		esc_data.zcd_intervals_pointer = 0;
}
