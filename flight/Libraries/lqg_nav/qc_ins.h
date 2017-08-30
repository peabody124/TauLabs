/**
 ******************************************************************************
 * @addtogroup TauLabsLibraries Tau Labs Libraries
 * @{
 * @addtogroup Quadcopter INS estimation
 * @{
 *
 * @file       qc_ins.h
 * @author     Tau Labs, http://taulabs.org, Copyright (C) 2017
 * @brief      Quadcopter INS estimation
 *
 * @see        The GNU Public License (GPL) Version 3
 *
 ******************************************************************************/
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

bool qcins_alloc(uintptr_t *qcins_handle);
bool qcins_init(uintptr_t qcins_handle);
bool qcins_predict(uintptr_t qcins_handle, const float roll, const float pitch, const float yaw, const float throttle, double Ts);
bool qcins_correct_accel_gyro(uintptr_t qcins_handle, const float gyros[3], const float accels[3]);
bool qcins_get_altitude(uintptr_t qcins_handle, float p[1]);
bool qcins_get_velocity(uintptr_t qcins_handle, float v[3]);
bool qcins_get_attitude(uintptr_t qcins_handle, float q[4]);
bool qcins_get_rate(uintptr_t qcins_handle, float rate[3]);
bool qcins_get_torque(uintptr_t qcins_handle, float torque[4]);

/**
 * @}
 * @}
 */