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
#include "esc_settings.h"

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

EscSettingsData config;

static void go_esc_nothing(uint16_t time) {};
// Fault and stopping transitions
static void go_esc_fault(uint16_t);
static void go_esc_stopping(uint16_t);
static void go_esc_stopped(uint16_t);

// Play sound on initially getting armed
static void go_esc_armed_sound(uint16_t);

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
			[ESC_EVENT_ARM] = ESC_STATE_ARMED_SOUND,
		},
	},
	[ESC_STATE_ARMED_SOUND] = {
		.entry_fn = go_esc_armed_sound,
		.next_state = {
			[ESC_EVENT_AUTO] = ESC_STATE_IDLE,
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

/**
 * Generate a tone by periodically applying power to one phase while
 * grounding two others
 *
 * @param duration how many ms to generate tone for
 * @param frequency frequency in hz to generate tone for
 */
static void esc_tone(uint32_t duration_ms, uint32_t frequency)
{
	const uint32_t on_duration_us = 10;
	const uint32_t period_us = 1e6 / frequency;
	
	// Want to be reasonable about duty cycle
	if(on_duration_us > period_us / 10)
		return;
	
	const uint32_t off_duration_us = period_us - on_duration_us;
	
	uint32_t start_time = PIOS_DELAY_GetRaw();
	PIOS_ESC_Arm();
	while(PIOS_DELAY_DiffuS(start_time) < duration_ms * 1000) {
		PIOS_ESC_BeepOn();
		PIOS_DELAY_WaituS(on_duration_us);
		
		PIOS_ESC_BeepOff();
		PIOS_DELAY_WaituS(off_duration_us);
		
		PIOS_WDG_UpdateFlag(1);
	}
	PIOS_ESC_Off();
}

void esc_process_static_fsm_rxn() {

	static uint32_t zero_time;
	static uint32_t stop_time;

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

			if(esc_data.current > config.SoftCurrentLimit) {
				esc_data.duty_cycle -= 1;
				if(esc_data.duty_cycle < config.MinDc)
					esc_data.duty_cycle = config.MinDc;

				PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
			}
			last_timer = cur_timer;

			if(esc_data.speed_setpoint == 0 && PIOS_DELAY_DiffuS(stop_time) > 100000) {
				esc_fsm_inject_event(ESC_EVENT_STOP,0);
			} else if (esc_data.speed_setpoint > 0)
				stop_time = PIOS_DELAY_GetRaw();
	
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

static void go_esc_armed_sound(uint16_t time)
{
	esc_tone(150, NOTES_AS);
	PIOS_DELAY_WaitmS(50);
	esc_tone(350, NOTES_C);
	PIOS_DELAY_WaitmS(50);
	esc_tone(600, NOTES_F0);
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
	esc_data.current_speed = config.InitialStartupSpeed;
	esc_data.consecutive_missed = 0;
	esc_data.consecutive_detected = 0;
	
	for(uint8_t i = 0; i < NUM_STORED_SWAP_INTERVALS; i++)
		esc_data.swap_intervals[i] = 0;
	esc_data.swap_interval_sum = 0;
	esc_data.duty_cycle = 0.12 * PIOS_ESC_MAX_DUTYCYCLE;
	PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, 50000);  // Grab stator for 50 ms
}

/**
 * Called each time we enter the startup wait state or when a commutation occurs after ZCD
 * schedule a commutation.
 */
static void go_esc_startup_wait(uint16_t time)
{
	commutate();
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed));
}

/**
 * ZCD detected. Hang out here till commutate.
 */
static void go_esc_startup_zcd(uint16_t time)
{
	if(PIOS_DELAY_DiffuS(esc_data.last_swap_time) > (RPM_TO_US(config.FinalStartupSpeed)/2) &&
	   PIOS_DELAY_DiffuS(esc_data.last_swap_time) < (RPM_TO_US(config.InitialStartupSpeed)/2)) {
		esc_data.consecutive_detected++;
		esc_data.consecutive_missed = 0;
		zcd(time);
	} else {
		esc_data.consecutive_detected = 0;
		esc_data.consecutive_missed++;
	}

	// Schedule next commutation but only if the last ZCD was near where we expected.  Essentially
	// this is dealing with a startup condition where noise will cause the ZCD to detect early, and
	// then drive the speed up and make this happen really quickly and trigger a stall.  It can be
	// neglected entirely but doing this makes it startup faster by phasing itself.
	if(abs(esc_data.last_zcd_time - RPM_TO_US(esc_data.current_speed) / 2) < 250)
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed) / 2);

	if(esc_data.consecutive_detected > 10) {
		esc_fsm_inject_event(ESC_EVENT_CLOSED, 0);
	} else {
		// TODO: Bring back current control loop
	}
}

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
		if(esc_data.current_speed < config.FinalStartupSpeed)
			esc_data.current_speed+=10;

		// Schedule next commutation
		esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, RPM_TO_US(esc_data.current_speed));
	} 
	
	// TODO: bring back current limit
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
uint32_t commutation_time;
static void go_esc_cl_commutated(uint16_t time)
{
	uint32_t timeval = PIOS_DELAY_GetRaw();
	commutate();
	esc_fsm_schedule_event(ESC_EVENT_TIMEOUT, esc_data.swap_interval_smoothed << 1);
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

	uint32_t zcd_delay = esc_data.swap_interval_smoothed * (30 - config.CommutationPhase) / 60 - config.CommutationOffset;
	esc_fsm_schedule_event(ESC_EVENT_COMMUTATED, zcd_delay);
	
	zcd1_time = PIOS_DELAY_DiffuS(timeval);
	
	if(esc_data.current > config.SoftCurrentLimit) {
		esc_data.duty_cycle -= config.MaxDcChange;
		bound_duty_cycle();
		PIOS_ESC_SetDutyCycle(esc_data.duty_cycle);
	} else {
#ifdef OPEN_LOOP
		esc_data.duty_cycle = esc_data.speed_setpoint * PIOS_ESC_MAX_DUTYCYCLE / 8000;
#else
		int32_t error = esc_data.speed_setpoint - esc_data.current_speed;

		esc_data.error_accum += error;
		if(esc_data.error_accum > config.ILim)
			esc_data.error_accum = config.ILim;
		if(esc_data.error_accum < -config.ILim)
			esc_data.error_accum = -config.ILim;

		zcd2_time = PIOS_DELAY_DiffuS(timeval);

		new_dc = (esc_data.speed_setpoint * config.Kff - config.Kff2 +
		      error * config.Kp + esc_data.error_accum * config.Ki) * PIOS_ESC_MAX_DUTYCYCLE / PID_SCALE; //* error + esc_data.error_accum;
		
		// For now keep this calculation as a float and rescale it here
		if((new_dc - esc_data.duty_cycle) > config.MaxDcChange)
			esc_data.duty_cycle += config.MaxDcChange;
		else if((new_dc - esc_data.duty_cycle) < -config.MaxDcChange)
			esc_data.duty_cycle -= config.MaxDcChange;
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

	if (esc_settings_load(&config) != 0) {
		esc_settings_defaults(&config);
	}
	
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

	esc_tone(300, NOTES_F);
	PIOS_DELAY_WaitmS(50);
	esc_tone(100, NOTES_F);
	PIOS_DELAY_WaitmS(50);	
	esc_tone(180, NOTES_C);
	PIOS_DELAY_WaitmS(20);
	esc_tone(1000, NOTES_F);

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
	//esc_data.zcd_fraction += ((float) esc_data.last_zcd_time / esc_data.last_swap_interval - esc_data.zcd_fraction) * 0.05;
	zcds++;
}

static void bound_duty_cycle()
{
	if(esc_data.duty_cycle < config.MinDc)
		esc_data.duty_cycle = config.MinDc;
	if(esc_data.duty_cycle > config.MaxDc)
		esc_data.duty_cycle = config.MaxDc;
}

