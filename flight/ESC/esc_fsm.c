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

// TOOD: Important safety features
// 1. Disable after a timeout from last control command
// 2. Reliable shutdown on high current
// 3. Monitor high voltage level to track battery (bonus for determining # of cells)
// 4. BEEPS!

static void go_esc_fault();
static void go_stopping();
static void go_stopped();
static void go_esc_init();
static void go_fsm_startup_grab();
static void go_fsm_startup_wait();
static void go_fsm_startup_zcd();
static void go_fsm_startup_nozcd();

static void esc_fsm_process_auto();
static void esc_fsm_inject_event();
static void esc_fsm_init();

enum esc_fsm_state {
	ESC_STATE_FSM_FAULT = 0,               /* Must be zero so undefined transitions land here */

	ESC_STATE_INIT,                        /* Perform self test, load settings */
	ESC_STATE_WAIT_FOR_ARM,                /* Cannot be used until throttle < threshold */

	ESC_STATE_IDLE,                        /* Waiting for throttle command */
	ESC_STATE_STOPPING,                    /* Safely shut down motor then wait for throttle */
	ESC_STATE_STOPPED,                     /* Motor shut down */

	ESC_STATE_STARTUP_GRAB,                /* Current limited hold state to align stator */
	ESC_STATE_STARTUP_WAIT,                /* Waiting for either a ZCD or commutation to occur */
	ESC_STATE_STARTUP_ZCD_DETECTED,        /* A ZCD detected, increment counter */
	ESC_STATE_STARTUP_NOZCD_COMMUTATED,    /* Commutation without a ZCD, reset counter */

	ESC_STATE_CL_COMMUTATED,
	ESC_STATE_CL_MISSED_ZCD,
	ESC_STATE_CL_ZCD,

	ESC_FSM_NUM_STATES
};

enum esc_event {
	ESC_EVENT_FAULT,
	ESC_EVENT_ARM,
	ESC_EVENT_START,
	ESC_EVENT_COMMUTATED,
	ESC_EVENT_ZCD,
	ESC_EVENT_CLOSED,
	ESC_EVENT_OVERCURRENT,
	ESC_EVENT_TIMEOUT,
	ESC_EVENT_AUTO,
	ESC_EVENT_NUM_EVENTS	/* Must be last */
};

struct esc_transition {
	void (*entry_fn) ();
	enum esc_state next_state[ESC_EVENT_NUM_EVENTS];
};

const static struct esc_transition esc_transition[ESC_FSM_NUM_STATES] = {
	[ESC_STATE_FSM_FAULT] = {
		.entry_fn = go_esc_fault,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_INIT] = {
		.entry_fn = go_esc_init,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_WAIT_FOR_ARM,
		},
	},
	[ESC_STATE_WAIT_FOR_ARM] = {
		.entry_fn = null,
		.next_state = {
			[ESC_EVENT_ARM] = ESC_STATE_IDLE,
		},
	},
	[ESC_STATE_IDLE] = {
		.entry_fn = null,
		.next_state = {
			[ESC_EVENT_START] = ESC_STATE_STARTUP_GRAB,
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
	[ESC_STATE_STARTUP_GRAB] = {
		.entry_fn = go_fsm_startup_grab,
		.next_state = {
			[ESC_EVENT_TIMEOUT] = ESC_STATE_STARTUP_WAIT,
		},
	},
	[ESC_STATE_STARTUP_WAIT] = {
		// Schedule commutation
		.entry_fn = go_fsm_startup_wait,
		.next_state = {
			// Could not get BEMF in time
			[ESC_EVENT_TIMEOUT] = ESC_STATE_STOPPING,
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_NOZCD_COMMUTATED,
		},
	},
	[ESC_STATE_STARTUP_ZCD_DETECTED] = {
		.entry_fn = go_fsm_startup_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_WAIT,
		},
	},
	[ESC_STATE_STARTUP_NOZCD_COMMUTATED] = {
		.entry_fn = go_fsm_startup_nozcd,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_STARTUP_WAIT,
		},
	},

	/* Closed loop states */
	[ESC_STATE_CL_COMMUTATED] = {
		.entry_fn = go_fsm_cl_commutated(),
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_CL_ZCD,
			[ESC_EVENT_TIMEOUT] = ESC_STATE_CL_MISSED_ZCD,
		},
	},
	[ESC_STATE_CL_MISSED_ZCD] = {
		.entry_fn = go_fsm_startup_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},
	[ESC_STATE_CL_ZCD] = {
		.entry_fn = go_fsm_startup_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
		},
	},

};

static enum esc_fsm_state esc_state;

void esc_process_static_fsm_rxn() {
	switch(esc_state) {
		case ESC_STATE_FSM_FAULT:
			// Shouldn't happen
			break;
		case ESC_STATE_INIT:
			// No static rxn?
			break;
		case ESC_STATE_WAIT_FOR_ARM:
			// Monitor for low input
			break;
		case ESC_STATE_IDLE:
			// Monitor for non-low input, start trigger?
			break;
		case ESC_STATE_STOPPING:
			// Monitor that the motor has stopped and current < threshold
			break;
		case ESC_STATE_STOPPED:
			// Not used
			break;
		case ESC_STATE_STARTUP_GRAB:
			// Small control loop for constant current
			break;
		case ESC_STATE_STARTUP_WAIT:
			// Timing adjusted in entry function
			// This should perform a current control (?)
			break;
		case ESC_STATE_STARTUP_ZCD_DETECTED:
			/* A ZCD detected, increment counter */
		case ESC_STATE_STARTUP_NOZCD_COMMUTATED:
			/* Commutation without a ZCD, reset counter */

		// TODO
		case ESC_STATE_CL_COMMUTATED:
		case ESC_STATE_CL_MISSED_ZCD:
		case ESC_STATE_CL_ZCD:
			break;
	}
}

/**
 * Initialize the FSM
 */
static void esc_fsm_init()
{
	pios_esc_off();
	esc_state = ESC_STATE_INIT;
}

/**
 * Invalid state transition occurred.
 */
static void go_esc_fault()
{
	pios_esc_off();
}

static void go_stopping()
{
	pios_esc_off();
}

static void go_stopped()
{
	pios_esc_off();
}

static void go_fsm_startup_grab()
{
}

/**
 * Called each time we enter the startup wait state
 * schedule a commutation and adjust time to cause the motor to accelerate
 * static reaction will handle fixed current
 */
static void go_fsm_startup_wait()
{
}


static void go_fsm_startup_zcd();
static void go_fsm_startup_nozcd();

static void esc_fsm_process_auto();
static void esc_fsm_inject_event();
