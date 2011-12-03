/**
 ******************************************************************************
 *
 * @file       esc.h
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      Main ESC header.
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

#ifndef ESC_H
#define ESC_H

/* PIOS Includes */
#include <pios.h>
#include "esc_fsm.h"
#include "esc_settings.h"
#include "esc_serial.h"
#include "esc_logger.h"

enum esc_control_method {
	ESC_CONTROL_PWM,
	ESC_CONTROL_SERIAL
};

enum esc_logging_status {
	ESC_LOGGING_NONE,
	ESC_LOGGING_CAPTURE,
	ESC_LOGGING_FULL,
	ESC_LOGGING_ECHO
};

struct esc_control {
	int16_t pwm_input;
	int16_t serial_input;
	enum esc_control_method control_method;
	bool serial_logging_enabled;
	enum esc_logging_status backbuffer_logging_status;
	bool save_requested;
};

#endif /* ESC_H */
