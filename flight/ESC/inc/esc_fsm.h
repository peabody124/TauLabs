/**
 ******************************************************************************
 *
 * @file       esc_fsm.h
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      Main ESC state machine header.
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

#ifndef ESC_FSM_H
#define ESC_FSM_H

#include "stdint.h"

#define NUM_STORED_SWAP_INTERVALS 6

enum esc_fsm_state {
	ESC_STATE_FSM_FAULT = 0,               /* Must be zero so undefined transitions land here */

	ESC_STATE_INIT,                        /* Perform self test, load settings */
	ESC_STATE_WAIT_FOR_ARM,                /* Cannot be used until throttle < threshold */
	ESC_STATE_ARMED_SOUND,                 /* Play the sound to indicate got throttle */
	
	ESC_STATE_IDLE,                        /* Waiting for throttle command */
	ESC_STATE_STOPPING,                    /* Safely shut down motor then wait for throttle */
	ESC_STATE_STOPPED,                     /* Motor shut down */

	ESC_STATE_STARTUP_ENABLE,              /* Start all the drivers and such */
	ESC_STATE_STARTUP_GRAB,                /* Current limited hold state to align stator */
	ESC_STATE_STARTUP_WAIT,                /* Waiting for either a ZCD or commutation to occur */
	ESC_STATE_STARTUP_ZCD_DETECTED,        /* A ZCD detected, increment counter */
	ESC_STATE_STARTUP_NOZCD_COMMUTATED,    /* Commutation without a ZCD, reset counter */

	ESC_STATE_CL_START,
	ESC_STATE_CL_COMMUTATED,
	ESC_STATE_CL_NOZCD,
	ESC_STATE_CL_ZCD,

	ESC_FSM_NUM_STATES
};

enum esc_event {
	ESC_EVENT_FAULT,
	ESC_EVENT_ARM,
	ESC_EVENT_START,
	ESC_EVENT_STOP,
	ESC_EVENT_COMMUTATED,
	ESC_EVENT_ZCD,
	ESC_EVENT_CLOSED,
	ESC_EVENT_OVERCURRENT,
	ESC_EVENT_LATE_COMMUTATION,
	ESC_EVENT_TIMEOUT,
	ESC_EVENT_AUTO,
	ESC_EVENT_NUM_EVENTS	/* Must be last */
};

struct esc_fsm_data {
	int32_t duty_cycle;
	int16_t current_speed;
	int16_t speed_setpoint;
	int32_t duty_cycle_setpoint;

	uint8_t scheduled_event_armed;
	uint16_t scheduled_event_time;
	enum esc_event scheduled_event;

	// Current and battery measurements
	uint32_t current_ma;
	uint32_t battery_mv;
	uint32_t smoothed_high_mv;
	
	uint32_t last_swap_time;
	uint32_t last_zcd_time;
	uint16_t last_swap_interval;
	uint16_t swap_intervals[NUM_STORED_SWAP_INTERVALS];
	uint32_t swap_interval_sum;
	uint16_t swap_interval_smoothed;
	uint8_t swap_intervals_pointer;
	uint16_t zcd_intervals[NUM_STORED_SWAP_INTERVALS];
	uint8_t zcd_intervals_pointer;
	int8_t zcd_fraction;
	uint16_t consecutive_detected;
	uint16_t consecutive_missed;
	uint16_t total_missed;
	uint16_t faults;
	uint8_t detected;
	int32_t error_accum;
	enum esc_fsm_state pre_fault_state;
	enum esc_event pre_fault_event;
	enum esc_fsm_state state;
	uint8_t locked;
};

#define PID_SCALE 32178

struct esc_transition {
	void (*entry_fn) (uint16_t time);
	enum esc_fsm_state next_state[ESC_EVENT_NUM_EVENTS];
};

struct esc_fsm_data * esc_fsm_init();
void esc_fsm_inject_event(enum esc_event event, uint16_t time);
void esc_process_static_fsm_rxn();


#define NOTES_F0 361
#define NOTES_FS 386
#define NOTES_G  406
#define NOTES_GS 422
#define NOTES_A  460
#define NOTES_AS 490
#define NOTES_B  520
#define NOTES_C  555
#define NOTES_CS 575
#define NOTES_D  618
#define NOTES_DS 665
#define NOTES_E  700
#define NOTES_F  745

#define ESC_CONFIG_MAGIC 0x763fedc

#endif
