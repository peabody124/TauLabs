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

// To act like a normal open loop ESC
//#define OPEN_LOOP

//#define RPM_TO_US(x) (1e6 * 60 / (x * COMMUTATIONS_PER_ROT) )
#define RPM_TO_US(x) (1428571 / x)
#define US_TO_RPM(x) RPM_TO_US(x)
#define COMMUTATIONS_PER_ROT (7*6)
#define ESC_CONFIG_MAGIC 0x763fedc

#define PID_SCALE 32178
struct esc_config config = {
	.max_dc_change = 0.1 * PIOS_ESC_MAX_DUTYCYCLE,
	.kp = 0.0001 * PID_SCALE,
	.ki = 0.00001 * PID_SCALE,
	.kff = 1.3e-4 * PID_SCALE,
	.kff2 = -0.05 * PID_SCALE,
	.ilim = 500,
	.min_dc = 0,
	.max_dc = 0.90 * PIOS_ESC_MAX_DUTYCYCLE,
	.initial_startup_speed = 100,
	.final_startup_speed = 700,
	.startup_current_target = 100,
	.commutation_phase = 17,
	.soft_current_limit = 250,
	.hard_current_limit = 600,
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
static void bound_duty_cycle();

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
			[ESC_EVENT_ARM] = ESC_STATE_IDLE,
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
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_GRAB,
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
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_STARTUP_ZCD_DETECTED] = {
		.entry_fn = go_esc_startup_zcd,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_WAIT,
			[ESC_EVENT_CLOSED] = ESC_STATE_CL_START,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_STARTUP_NOZCD_COMMUTATED] = {
		.entry_fn = go_esc_startup_nozcd,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_STARTUP_ZCD_DETECTED,
			[ESC_EVENT_COMMUTATED] = ESC_STATE_STARTUP_NOZCD_COMMUTATED,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},

	/* Closed loop states */
	[ESC_STATE_CL_START] = {
		.entry_fn = go_esc_cl_start,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_CL_COMMUTATED] = {
		.entry_fn = go_esc_cl_commutated,
		.next_state = {
			[ESC_EVENT_ZCD] = ESC_STATE_CL_ZCD,
			[ESC_EVENT_TIMEOUT] = ESC_STATE_CL_NOZCD,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_CL_NOZCD] = {
		.entry_fn = go_esc_cl_nozcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
	[ESC_STATE_CL_ZCD] = {
		.entry_fn = go_esc_cl_zcd,
		.next_state = {
			[ESC_EVENT_COMMUTATED] = ESC_STATE_CL_COMMUTATED,
			[ESC_EVENT_LATE_COMMUTATION] = ESC_STATE_CL_COMMUTATED,
			[ESC_EVENT_STOP] = ESC_STATE_STOPPING,
		},
	},
};

static struct esc_fsm_data esc_data;
uint32_t stops_requested, stops_from_here;
void esc_process_static_fsm_rxn() {

	static uint32_t zero_time;

	switch(esc_data.state) {
		case ESC_STATE_FSM_FAULT:
			// Shouldn't happen
			break;
		case ESC_STATE_INIT:
			// No static rxn.  Goes straight to wait_for_arm.
			break;
		case ESC_STATE_WAIT_FOR_ARM:
		case ESC_STATE_STOPPED:
			if(esc_data.speed_setpoint == 0 && zero_time > 10000) {
				esc_fsm_inject_event(ESC_EVENT_ARM,0);
				zero_time = 0;
			}
			else zero_time++;
			break;
		case ESC_STATE_IDLE:
			if(esc_data.speed_setpoint > 0)
				esc_fsm_inject_event(ESC_EVENT_START,0);
			break;
		case ESC_STATE_STOPPING:
			// Monitor that the motor has stopped and current < threshold
			PIOS_ESC_Off();
			break;

		/**
		 * Startup state transitions
		 */
		case ESC_STATE_STARTUP_GRAB:
		case ESC_STATE_STARTUP_WAIT:
		case ESC_STATE_STARTUP_ZCD_DETECTED:
		case ESC_STATE_STARTUP_NOZCD_COMMUTATED:			
			if(esc_data.speed_setpoint == 0) {
				stops_from_here++;
				esc_fsm_inject_event(ESC_EVENT_STOP,0);
			}
			
			break;

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
				esc_data.duty_cycle -= 1;
				if(esc_data.duty_cycle < config.min_dc)
					esc_data.duty_cycle = config.min_dc;

				PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
			}
			last_timer = cur_timer;

			if(esc_data.speed_setpoint == 0) {
				stops_requested++;
				esc_fsm_inject_event(ESC_EVENT_STOP,0);
			}
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
//	PIOS_ADC_StopDma();
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
	esc_data.consecutive_detected = 0;
	esc_data.current_speed = 0;
	esc_data.detected = true;
	esc_data.scheduled_event_armed = false;
}

/**
 * Enable ADC, timers, drivers, etc
 */
static void go_esc_startup_enable(uint16_t time)
{
	PIOS_ADC_StartDma();
	PIOS_ESC_SetMode(ESC_MODE_LOW_ON_PWM_HIGH);

	PIOS_ESC_Arm();
}

/**
 * Grab the stator for startup sequencing
 */
static void go_esc_startup_grab(uint16_t time)
{
	// TODO: Set up a timeout for whole startup system

	PIOS_ESC_SetState(0);
	esc_data.current_speed = config.initial_startup_speed;
	esc_data.consecutive_missed = 0;
	esc_data.consecutive_detected = 0;
	
	for(uint8_t i = 0; i < NUM_STORED_SWAP_INTERVALS; i++)
		esc_data.swap_intervals[i] = 0;
	esc_data.swap_interval_sum = 0;
	esc_data.duty_cycle = 0.10 * PIOS_ESC_MAX_DUTYCYCLE;
	PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, 30000);  // Grab stator for 30 ms
}

/**
 * Called each time we enter the startup wait state or when a commutation occurs after ZCD
 * schedule a commutation.
 */
uint32_t startup_waits = 0;
static void go_esc_startup_wait(uint16_t time)
{
	commutate();
	startup_waits ++;
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed));
}

/**
 * ZCD detected. Hang out here till commutate.
 */
 
uint16_t current_time;
uint16_t prev_time;
int16_t diff_time;
static void go_esc_startup_zcd(uint16_t time)
{

	esc_data.consecutive_detected++;
	esc_data.consecutive_missed = 0;

	zcd(time);

	// Since we aren't getting ZCD keep accelerating
	if(esc_data.current_speed < config.final_startup_speed)
		esc_data.current_speed+=5;

	// Schedule next commutation but only if the last ZCD was near where we expected.  Essentially
	// this is dealing with a startup condition where noise will cause the ZCD to detect early, and
	// then drive the speed up and make this happen really quickly and trigger a stall.  It can be
	// neglected entirely but doing this makes it startup faster by phasing itself.
	if(abs(esc_data.last_zcd_time - RPM_TO_US(esc_data.current_speed)) < 250)
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed) / 2);

	prev_time = TIM4->CCR1;
	current_time = TIM4->CNT + RPM_TO_US(esc_data.current_speed) / 2;
	diff_time = prev_time - current_time;
	
	if(esc_data.consecutive_detected > 6) {
		esc_fsm_inject_event(ESC_EVENT_CLOSED, 0);
	} else {
		// Timing adjusted in entry function
		// This should perform run a current control loop

/*		float current_error = (config.startup_current_target - esc_data.current);
		current_error *= 0.00001;
		esc_data.duty_cycle += current_error;

		if(esc_data.duty_cycle < 0.01)
			esc_data.duty_cycle = 0.01;
		if(esc_data.duty_cycle > 0.4)
			esc_data.duty_cycle = 0.4;
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle); */
	}
}

uint32_t startup_schedules = 0;
/**
 * Commutation occured during startup without a ZCD
 */
static void go_esc_startup_nozcd(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed++;

	if(esc_data.consecutive_missed > 100) {
		esc_fsm_inject_event(ESC_EVENT_FAULT, time);
	} else {
		commutate();

		// Since we aren't getting ZCD keep accelerating
		if(esc_data.current_speed < config.final_startup_speed)
			esc_data.current_speed+=5;

		startup_schedules++;
		
		// Schedule next commutation
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed));
	} 
			
	// Timing adjusted in entry function
	// This should perform run a current control loop
	int32_t current_error = (config.startup_current_target - esc_data.current);
	current_error /= 8;
	esc_data.duty_cycle += current_error;

	bound_duty_cycle();
	PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
}

uint32_t starts = 0;

/**
 * Commutation in closed loop
 */
static void go_esc_cl_start(uint16_t time)
{
	starts++;
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed = 0;
}

/**
 * Commutation in closed loop
 */
uint32_t commutation_time;
static void go_esc_cl_commutated(uint16_t time)
{
	uint32_t timeval = PIOS_DELAY_GetRaw();
	commutate();
	esc_fsm_schedule_event(ESC_EVENT_TIMEOUT, esc_data.swap_interval_smoothed * 1.5);
//	esc_data.Kv += (esc_data.current_speed / (12 * esc_data.duty_cycle) - esc_data.Kv) * 0.001;
	commutation_time = PIOS_DELAY_DiffuS(timeval);
}

/**
 * When a zcd is detected
 */
uint32_t zcd1_time;
uint32_t zcd2_time;
uint32_t zcd3_time;
uint32_t zcd4_time;
uint32_t zcd5_time;
int32_t new_dc;
static void go_esc_cl_zcd(uint16_t time)
{
	uint32_t timeval = PIOS_DELAY_GetRaw();
	
	esc_data.consecutive_detected++;
	esc_data.consecutive_missed = 0;

	zcd(time);

	uint32_t zcd_delay = esc_data.swap_interval_smoothed * (30 - config.commutation_phase) / 60;
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, zcd_delay);
	
	zcd1_time = PIOS_DELAY_DiffuS(timeval);
	
	if(esc_data.current > config.soft_current_limit) {
		esc_data.duty_cycle -= config.max_dc_change;
		bound_duty_cycle();
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	} else {
#ifdef OPEN_LOOP
		esc_data.duty_cycle = esc_data.speed_setpoint * PIOS_ESC_MAX_DUTYCYCLE / 8000;
#else
		int32_t error = esc_data.speed_setpoint - esc_data.current_speed;

		esc_data.error_accum += error;
		if(esc_data.error_accum > config.ilim)
			esc_data.error_accum = config.ilim;
		if(esc_data.error_accum < -config.ilim)
			esc_data.error_accum = -config.ilim;

		zcd2_time = PIOS_DELAY_DiffuS(timeval);

		new_dc = (esc_data.speed_setpoint * config.kff - config.kff2 +
		      error * config.kp + esc_data.error_accum * config.ki) * PIOS_ESC_MAX_DUTYCYCLE / PID_SCALE; //* error + esc_data.error_accum;
		
		// For now keep this calculation as a float and rescale it here
		if((new_dc - esc_data.duty_cycle) > config.max_dc_change)
			esc_data.duty_cycle += config.max_dc_change;
		else if((new_dc - esc_data.duty_cycle) < -config.max_dc_change)
			esc_data.duty_cycle -= config.max_dc_change;
		else
			esc_data.duty_cycle = new_dc;

		zcd3_time = PIOS_DELAY_DiffuS(timeval);

#endif	
		bound_duty_cycle();
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	}
	zcd4_time = PIOS_DELAY_DiffuS(timeval);
}

/**
 * For now fault out if we miss a zcd
 */
static void go_esc_cl_nozcd(uint16_t time)
{
	esc_data.consecutive_detected = 0;
	esc_data.consecutive_missed++;
	if(esc_data.consecutive_missed > 10000)
		esc_fsm_inject_event(ESC_EVENT_FAULT, 0);
	else {
//		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle / 10);
		esc_fsm_inject_event(ESC_EVENT_COMMUTATED, 0);
	}
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
	esc_data.detected = true;

	esc_fsm_process_auto(PIOS_DELAY_GetuS());

	return &esc_data;
}


/**
 * Respond to an event with state system
 */
enum esc_event last_event;

void esc_fsm_inject_event(enum esc_event event, uint16_t time)
{
	PIOS_IRQ_Disable();

	if(esc_transition[esc_data.state].next_state[event] == ESC_STATE_FSM_FAULT) {
		esc_data.pre_fault_state = esc_data.state;
		esc_data.pre_fault_event = event;
	}

	/*
	 * Move to the next state
	 *
	 * This is done prior to calling the new state's entry function to
	 * guarantee that the entry function never depends on the previous
	 * state.  This way, it cannot ever know what the previous state was.
	 */
	esc_data.state = esc_transition[esc_data.state].next_state[event];
	last_event = event;
	
	PIOS_IRQ_Enable();

	/* Call the entry function (if any) for the next state. */
	if (esc_transition[esc_data.state].entry_fn) {
		esc_transition[esc_data.state].entry_fn(time);
	}

	/* Process any AUTO transitions in the FSM */
	esc_fsm_process_auto(time);
}

/**
 * Process any automatic transitions
 */
static void esc_fsm_process_auto(uint16_t time)
{

	while (esc_transition[esc_data.state].next_state[ESC_EVENT_AUTO]) {
		PIOS_IRQ_Disable();
		esc_data.state = esc_transition[esc_data.state].next_state[ESC_EVENT_AUTO];
		PIOS_IRQ_Enable();

		/* Call the entry function (if any) for the next state. */
		if (esc_transition[esc_data.state].entry_fn) {
			esc_transition[esc_data.state].entry_fn(time);
		}
	}

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
	esc_data.scheduled_event_time = time + TIM4->CNT;
	PIOS_IRQ_Enable();
	TIM_SetCompare1(TIM4, esc_data.scheduled_event_time);
}

uint32_t timeouts = 0;
uint32_t events = 0;
/**
 * The event interrupt
 */
void PIOS_DELAY_timeout() {
	timeouts++;
	if(esc_data.scheduled_event_armed) {
		events++;
		esc_data.scheduled_event_armed = false;
		esc_fsm_inject_event(esc_data.scheduled_event, PIOS_DELAY_GetuS());
	}
}

uint32_t commutations;
/**
 * Function to commutate the ESC and store any necessary information
 */
static void commutate()
{
	esc_data.last_swap_interval = PIOS_DELAY_DiffuS(esc_data.last_swap_time);
	esc_data.last_swap_time = PIOS_DELAY_GetRaw();

	esc_data.swap_interval_sum += esc_data.last_swap_interval - esc_data.swap_intervals[esc_data.swap_intervals_pointer];
	esc_data.swap_intervals[esc_data.swap_intervals_pointer] = esc_data.last_swap_interval;
	esc_data.swap_intervals_pointer++;
	if(esc_data.swap_intervals_pointer >= NUM_STORED_SWAP_INTERVALS)
		esc_data.swap_intervals_pointer = 0;	
	
	esc_data.swap_interval_smoothed = esc_data.swap_interval_sum / NUM_STORED_SWAP_INTERVALS;
	esc_data.current_speed = US_TO_RPM(esc_data.swap_interval_smoothed);
	
	esc_data.detected = false;
	commutations++;
	
	PIOS_ESC_NextState();
}

uint32_t zcds;
/**
 * Function to handle the book-keeping when a zcd detected
 */
static void zcd(uint16_t time)
{
	esc_data.last_zcd_time = PIOS_DELAY_DiffuS(esc_data.last_swap_time);
	esc_data.detected = true;
	esc_data.zcd_fraction += ((float) esc_data.last_zcd_time / esc_data.last_swap_interval - esc_data.zcd_fraction) * 0.05;
	zcds++;
}

static void bound_duty_cycle()
{
	if(esc_data.duty_cycle < config.min_dc)
		esc_data.duty_cycle = config.min_dc;
	if(esc_data.duty_cycle > config.max_dc)
		esc_data.duty_cycle = config.max_dc;
}

