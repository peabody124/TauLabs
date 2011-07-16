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
	ESC_EVENT_COMMUTATED,
	ESC_EVENT_ZCD,
	ESC_EVENT_CLOSED,
	ESC_EVENT_OVERCURRENT,
	ESC_EVENT_TIMEOUT,
	ESC_EVENT_AUTO,
	ESC_EVENT_NUM_EVENTS	/* Must be last */
};

struct esc_fsm_data {
	float duty_cycle;
	float current;
	uint16_t current_speed;
	uint16_t speed_setpoint;
	uint8_t scheduled_event_armed;
	uint16_t last_swap_time;
	uint16_t last_zcd_time;
	uint16_t swap_intervals[NUM_STORED_SWAP_INTERVALS];
	uint8_t swap_intervals_pointer;
	uint16_t zcd_intervals[NUM_STORED_SWAP_INTERVALS];
	uint8_t zcd_intervals_pointer;
	uint16_t consecutive_detected;
	uint16_t consecutive_missed;
	uint16_t faults;
	uint8_t detected;
	float error_accum;
	float Kv;
	enum esc_fsm_state state;
	enum esc_event scheduled_event;
};

struct esc_config {
	float kp;
	float ki;
	float kff;
	float kff2;
	float accum;
	float ilim;
	float max_dc_change;
	float min_dc;
	float max_dc;
	uint16_t initial_startup_speed;
	uint16_t final_startup_speed;
	float commutation_phase;
	float soft_current_limit;
	float hard_current_limit;
	uint32_t magic;
};


struct esc_transition {
	void (*entry_fn) (uint16_t time);
	enum esc_fsm_state next_state[ESC_EVENT_NUM_EVENTS];
};

struct esc_fsm_data * esc_fsm_init();
void esc_fsm_inject_event(enum esc_event event, uint16_t time);
void esc_process_static_fsm_rxn();

#endif
