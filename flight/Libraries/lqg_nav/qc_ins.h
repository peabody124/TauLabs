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

bool qcins_predict(uintptr_t qcins_handle, const float roll, const float pitch, const float yaw, const float throttle, float Ts);

bool qcins_correct_accel_gyro(uintptr_t qcins_handle, const float accels[3], const float gyros[3]);
bool qcins_correct_accel_gyro(uintptr_t qcins_handle, const float accels[3], const float gyros[3]);
bool qcins_correct_baro(uintptr_t qcins_handle, float baro);
bool qcins_correct_mag(uintptr_t qcins_handle, const float mag[3]);

// Ideally these should be good defaults and not need adjusting
bool qcins_set_sensor_noise(uintptr_t qcins_handle, const float noises[9]);
bool qcins_set_process_noise(uintptr_t qcins_handle, const float noises[15]);

bool qcins_set_gains(uintptr_t qcins_handle, const float gains_new[4]);
bool qcins_set_tau(uintptr_t qcins_handle, const float tau_new);
bool qcins_set_mu(uintptr_t qcins_handle, const float mu_new);

bool qcins_get_altitude(uintptr_t qcins_handle, float p[1]);
bool qcins_get_velocity(uintptr_t qcins_handle, float v[3]);
bool qcins_get_attitude(uintptr_t qcins_handle, float q[4]);
bool qcins_get_rate(uintptr_t qcins_handle, float rate[3]);
bool qcins_get_torque(uintptr_t qcins_handle, float torque[4]);

/**
 * @}
 * @}
 */
