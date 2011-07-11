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

#define NUM_STORED_SWAP_INTERVALS 6

#define RPM_TO_US(x) (1e6 * 60 / (x * COMMUTATIONS_PER_ROT) )
#define US_TO_RPM(x) RPM_TO_US(x)
#define COMMUTATIONS_PER_ROT (7*6)

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
static void go_esc_cl_commutated(uint16_t time) {};

// Misc utility functions
static void esc_fsm_process_auto(uint16_t);
static void esc_fsm_schedule_event(enum esc_event event, uint16_t time);

static void commutate();

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
			[ESC_EVENT_AUTO] = ESC_STATE_IDLE,
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
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_WAIT,
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
	[ESC_STATE_CL_COMMUTATED] = {
		.entry_fn = go_esc_cl_commutated,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_CL_ZCD,
			[ESC_EVENT_TIMEOUT] = ESC_STATE_CL_MISSED_ZCD,
		},
	},
	[ESC_STATE_CL_MISSED_ZCD] = {
		.entry_fn = go_esc_startup_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},
	[ESC_STATE_CL_ZCD] = {
		.entry_fn = go_esc_startup_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},

};

struct esc_fsm_data {
	float duty_cycle;
	uint8_t scheduled_event_armed;
	uint16_t last_swap_time;
	uint16_t swap_intervals[NUM_STORED_SWAP_INTERVALS];
	uint8_t swap_intervals_pointer;
	uint16_t consecutive_detected;
	uint16_t consecutive_missed;
	uint16_t current_speed;
	uint16_t faults;
	enum esc_fsm_state state;
	enum esc_event scheduled_event;
};

static struct esc_fsm_data esc_data;

void esc_process_static_fsm_rxn() {

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
			// Timing adjusted in entry function
			// This should perform run a current control loop
			break;

		/**
		 * Closed loop transitions
		 */
		case ESC_STATE_CL_COMMUTATED:
		case ESC_STATE_CL_MISSED_ZCD:
		case ESC_STATE_CL_ZCD:
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
	PIOS_ESC_Off();
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
	esc_data.current_speed = 100;
	esc_data.duty_cycle = 0.14;
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
/*	if(esc_data.consecutive_detected > THRESHOLD)
	inject event to go to closed loop
 */
}

/**
 * ZCD detected. Hang out here till commutate.
 */
static void go_esc_startup_zcd(uint16_t time)
{
	esc_data.consecutive_detected++;
	esc_data.consecutive_missed = 0;

	// Schedule next commutation
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, PIOS_DELAY_GetuS() + RPM_TO_US(esc_data.current_speed));
}

/**
 * Commutation occured during startup without a ZCD
 */
static void go_esc_startup_nozcd(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed++;
	commutate();

	// Since we aren't getting ZCD keep accelerating
	if(esc_data.current_speed < 2000)
		esc_data.current_speed+=2;

	// Schedule next commutation
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, PIOS_DELAY_GetuS() + RPM_TO_US(esc_data.current_speed));
}


/**
 * Initialize the FSM
 */
void esc_fsm_init()
{
	PIOS_ESC_Off();
	esc_data.state = ESC_STATE_INIT;
	esc_data.duty_cycle = 0;
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed = 0;
	esc_data.current_speed = 0;
	for (uint8_t i = 0; i < NUM_STORED_SWAP_INTERVALS; i++)
		esc_data.swap_intervals[i] = 0;
	esc_data.swap_intervals_pointer = 0;
	esc_data.last_swap_time = 0;
	esc_data.scheduled_event_armed = false;
	esc_data.faults = 0;

	esc_fsm_process_auto(PIOS_DELAY_GetuS());
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
	uint16_t swap_diff = (uint16_t) (this_swap_time - esc_data.last_swap_time);
	esc_data.last_swap_time = this_swap_time;

	esc_data.last_swap_time = this_swap_time;
	esc_data.swap_intervals[esc_data.swap_intervals_pointer++] = swap_diff;
	if(esc_data.swap_intervals_pointer > NUM_STORED_SWAP_INTERVALS)
		esc_data.swap_intervals_pointer = 0;

	PIOS_ESC_NextState();
}
