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

#include "pios_heap.h"
#include "math.h"
#include "physical_constants.h"
#include "stdint.h"

#include "stdio.h"

// Local definitions
#define NUMX 19  // number of state variables
#define NUMF 63  // number of non-zero variables stored in the F matrix
#define NUMH 32  // number of non-zero variables stored in the H matrix
#define NUMP 97  // number of non-zero variables in the upper triangular of covariance
#define NUMO 10  // number of total outputs

// Private methods
void covariance_init(float *Pnew);
void state_prediction(const float * restrict x, const float *u, float Ts, const float *restrict param, float * restrict xnew);
void covariance_prediction(const float *restrict P, const float * restrict F, const float * restrict Q, float * restrict Pnew);
void linearize_FH(const float * restrict x, float Ts, const float * restrict param, float * restrict F, float * restrict H);
void baro_correction(const float *restrict x, const float *restrict P, float baro,
        const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew);
void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, const float * restrict param,
        float *restrict xnew, float *restrict Pnew);
void gyro_correction(const float * restrict x, const float * restrict P, const float * restrict gyro,
	const float *restrict H, const float *restrict R,	float *restrict xnew, float *restrict Pnew);
void accel_correction(const float *restrict x, const float *restrict P, const float *restrict accel,
        const float *restrict H, const float *restrict R, const float *restrict param, float *restrict xnew, float *restrict Pnew);
void update_state(const float * restrict x, const float *restrict P, float *restrict xnew, float *restrict Pnew);
void update_state_z(const float * restrict xold, float *restrict xnew);
void normalize_state(float *x);

enum qcins_state_magic {
  QCINS_STATE_MAGIC = 0xe7b5500e, // echo "qc_ins.c" | md5
};

struct qcins_state {

	float Ts;      // The kalman filter time step

	struct {
		float g;
		float beta_r;
		float beta_p;
		float beta_y1;
		float beta_y2;
		float tau;
		float mu;
		float Be[3];
	} __attribute__((__packed__)) params;

	struct {
		float beta_t;
		float bias[3];
	} __attribute__((__packed__)) init;

	bool armed;

	float x[NUMX]; // buffer to store the current state
	float P[NUMP]; // buffer to store the covariance
	float F[NUMF]; // buffer to store the state dynamics derivatives
	float H[NUMH]; // buffer to store the output derivatives
	float Q[NUMX]; // process covariance
	float R[NUMO]; // measurement noises

	float xnew[NUMX];
	float Pnew[NUMP];

	enum qcins_state_magic magic;

};

bool qcins_validate(struct qcins_state *qcins_state)
{
	if (qcins_state == NULL)
		return false;

	return (qcins_state->magic == QCINS_STATE_MAGIC);
}

bool qcins_init(uintptr_t qcins_handle)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	// Initialize state and covariance
	for (uint32_t i = 0; i < NUMX; i++)
		qcins_state->x[i] = 0;
	qcins_state->x[4] = 1.0f;
	qcins_state->x[14] = 1.0f; // Initial thrust 1
	qcins_state->x[18] = 2.0f; // Initial thrust to weight ratio 2:1

	covariance_init(qcins_state->P);

	// Store the process noises. Very low except for drift in torques
	float *Q = qcins_state->Q;
	for (uint32_t i = 0; i < NUMX; i++)
		Q[i] = 1e-4;
	Q[0] = 1e-3;                        // Position
	Q[1] = Q[2] = 1e-3;                 // Horizontal velocity
	Q[3] = 1e-4;                        // Vertical velocity
	Q[4] = Q[5] = Q[6] = Q[7] = 1e-3;   // Quaternion
	Q[8] = Q[9] = Q[10] = 1e-1;         // Rotation rates
	Q[11] = Q[12] = Q[13] = 1;          // Torques
	Q[14] = 1e-1;                       // Thrust
	Q[15] = Q[16] = Q[17] = 1e-5f;      // Bias states
	Q[18] = 1e-5f;                      // Thrust gain

	qcins_state->armed = false;

	// Store the observation noises
	qcins_state->R[0] = 100;
	qcins_state->R[1] = qcins_state->R[2] = qcins_state->R[3] = 1e6;
	qcins_state->R[4] = qcins_state->R[5] = qcins_state->R[6] = 1e3;
	qcins_state->R[7] = qcins_state->R[8] = qcins_state->R[9] = 1e3;

	qcins_state->init.bias[0] = 0.0f;
	qcins_state->init.bias[1] = 0.0f;
	qcins_state->init.bias[2] = 0.0f;
	qcins_state->init.beta_t = 2.0f; // thrust to weight ratio 2:1

	// Defaults for the parameters
	qcins_state->params.g = 9.81f;
	qcins_state->params.beta_r = 10000.0f * DEG2RAD;
	qcins_state->params.beta_p = 10000.0f * DEG2RAD;
	qcins_state->params.beta_y1 = 1000.0f * DEG2RAD;
	qcins_state->params.beta_y2 = 1000.0f * DEG2RAD;
	qcins_state->params.tau = 0.050f;
	qcins_state->params.mu = 1;
	qcins_state->params.Be[0] = 1;
	qcins_state->params.Be[1] = 0;
	qcins_state->params.Be[2] = 0;
	
	return true;
}

bool qcins_alloc(uintptr_t *qcins_handle)
{
	struct qcins_state *qcins_state = (struct qcins_state *) PIOS_malloc(sizeof(struct qcins_state));
	if (qcins_state == NULL)
		return false;

	qcins_state->magic = QCINS_STATE_MAGIC;

	qcins_init((uintptr_t) qcins_state);

	(*qcins_handle) = (uintptr_t) qcins_state;

	return true;
}

// Ideally these should be good defaults and not need adjusting
bool qcins_set_sensor_noise(uintptr_t qcins_handle, const float noises[9])
{
   struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
   if (!qcins_validate(qcins_state))
      return false;

   for (int i = 0; i < NUMO; i++)
      qcins_state->R[i] = noises[i];
   return true;
}

// Ideally these should be good defaults and not need adjusting
bool qcins_set_process_noise(uintptr_t qcins_handle, const float noises[15])
{
   struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
   if (!qcins_validate(qcins_state))
      return false;

   for (int i = 0; i < NUMX; i++)
      qcins_state->Q[i] = noises[i];
   return true;
}

bool qcins_set_init_thrust(uintptr_t qcins_handle, const float beta_t_new)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->init.beta_t = beta_t_new;
	return true;
}

bool qcins_set_init_bias(uintptr_t qcins_handle, const float bias_new[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->init.bias[0] = bias_new[0];
	qcins_state->init.bias[1] = bias_new[1];
	qcins_state->init.bias[2] = bias_new[2];
	return true;
}

bool qcins_set_armed(uintptr_t qcins_handle, bool armed)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->armed = armed;
	return true;
}

bool qcins_set_gains(uintptr_t qcins_handle, const float gains_new[4])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->params.beta_r = expf(gains_new[0]) * DEG2RAD;
	qcins_state->params.beta_p = expf(gains_new[1]) * DEG2RAD;
	qcins_state->params.beta_y1 = expf(gains_new[2]) * DEG2RAD;
	qcins_state->params.beta_y2 = expf(gains_new[3]) * DEG2RAD;

	return true;
}

bool qcins_set_tau(uintptr_t qcins_handle, const float tau_new)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	// Expects to use the tau from system ident which is in log representation
	qcins_state->params.tau = expf(tau_new);
	return true;
}

bool qcins_set_mu(uintptr_t qcins_handle, const float mu_new)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->params.mu = mu_new;
	return true;
}

static void qcins_unarmed_state(struct qcins_state *qcins_state)
{
	// Set the torques to zero
	qcins_state->x[11] = 0.0f;
	qcins_state->x[12] = 0.0f;
	qcins_state->x[13] = 0.0f;
	// qcins_state->x[14] = 1.0f;

	// Set these system parameters to their initial values
	qcins_state->x[15] = qcins_state->init.bias[0];
	qcins_state->x[16] = qcins_state->init.bias[1];
	qcins_state->x[17] = qcins_state->init.bias[2];
	qcins_state->x[18] = qcins_state->init.beta_t; // thrust to weight ratio
}

bool qcins_predict(uintptr_t qcins_handle, const float roll, const float pitch, const float yaw, const float throttle, float Ts)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	/*printf("Params: %f %f %f %f %f %f %f %f %f %f %f\r\n",
		qcins_state->params.g, qcins_state->params.beta_r,
		qcins_state->params.beta_p, qcins_state->params.beta_y1,
		qcins_state->params.beta_y2, qcins_state->x[18],
		qcins_state->params.tau, qcins_state->params.mu,
		qcins_state->params.Be[0], qcins_state->params.Be[1], 
		qcins_state->params.Be[2]);*/

	const float u[4] = {roll, pitch, yaw, throttle};

	if (!qcins_state->armed) {
		qcins_unarmed_state(qcins_state);
	}

	// Take the state and covariance forward one time step
	linearize_FH(qcins_state->x, Ts, &qcins_state->params.g, qcins_state->F, qcins_state->H);
	state_prediction(qcins_state->x, u, Ts, &qcins_state->params.g, qcins_state->xnew);
	covariance_prediction(qcins_state->P, qcins_state->F, qcins_state->Q, qcins_state->Pnew);
	update_state(qcins_state->xnew, qcins_state->Pnew, qcins_state->x, qcins_state->P);
	normalize_state(qcins_state->x);

	return true;
}

bool qcins_correct_accel_gyro(uintptr_t qcins_handle, const float accels[3], const float gyros[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	gyro_correction(qcins_state->x, qcins_state->P, gyros, qcins_state->H, qcins_state->R, qcins_state->xnew, qcins_state->Pnew);
	accel_correction(qcins_state->xnew, qcins_state->Pnew, accels, qcins_state->H, qcins_state->R, &qcins_state->params.g, qcins_state->x, qcins_state->P);
	normalize_state(qcins_state->x);
	return true;
}

bool qcins_correct_baro(uintptr_t qcins_handle, float baro)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	baro_correction(qcins_state->x, qcins_state->P, baro, qcins_state->H, qcins_state->R, qcins_state->xnew, qcins_state->Pnew);
	update_state(qcins_state->xnew, qcins_state->Pnew, qcins_state->x, qcins_state->P);
	normalize_state(qcins_state->x);
	return true;
}

bool qcins_correct_mag(uintptr_t qcins_handle, const float mag[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	// Only preserve direction
	const float mag_len = sqrtf(mag[0]*mag[0] + mag[1]*mag[1]) + 1e-10f;
	const float mag_norm[2] = {mag[0] / mag_len, mag[1] / mag_len};

	mag_correction(qcins_state->x, qcins_state->P, mag_norm, qcins_state->H, qcins_state->R, &qcins_state->params.g, qcins_state->xnew, qcins_state->Pnew);
	update_state_z(qcins_state->x, qcins_state->xnew); // Only keep change in yaw
	update_state(qcins_state->xnew, qcins_state->Pnew, qcins_state->x, qcins_state->P);
	normalize_state(qcins_state->x);
	return true;
}

/**
 * @param[in] qcins_handle handle for estimation
 * @param[out] p store the current altitude here
 */
bool qcins_get_altitude(uintptr_t qcins_handle, float *p)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	p[0] = qcins_state->x[0];
	return true;
}

/**
 * @param[in] qcins_handle handle for estimation
 * @param[out] v store the current velocity here
 */
bool qcins_get_velocity(uintptr_t qcins_handle, float v[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	v[0] = qcins_state->x[1];
	v[1] = qcins_state->x[2];
	v[2] = qcins_state->x[3];
	return true;
}

/**
 * @param[in] qcins_handle handle for estimation
 * @param[out] q store the current attitude here
 */
bool qcins_get_attitude(uintptr_t qcins_handle, float q[4])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	q[0] = qcins_state->x[4];
	q[1] = qcins_state->x[5];
	q[2] = qcins_state->x[6];
	q[3] = qcins_state->x[7];

	return true;
}

/**
 * @param[in] qcins_handle handle for estimation
 * @param[out] rate store the current rate here
 */
bool qcins_get_rate(uintptr_t qcins_handle, float rate[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	rate[0] = qcins_state->x[8] * 180.0f / (float) M_PI;
	rate[1] = qcins_state->x[9] * 180.0f / (float) M_PI;
	rate[2] = qcins_state->x[10] * 180.0f / (float) M_PI;
	return true;
}

/**
 * @param[in] qcins_handle handle for estimation
 * @param[out] torque store the current torques here including thrust
 */
bool qcins_get_torque(uintptr_t qcins_handle, float torque[4])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	torque[0] = qcins_state->x[11];
	torque[1] = qcins_state->x[12];
	torque[2] = qcins_state->x[13];
	torque[3] = qcins_state->x[14];
	return true;
}

bool qcins_get_bias(uintptr_t qcins_handle, float bias[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	bias[0] = qcins_state->x[15];
	bias[1] = qcins_state->x[16];
	bias[2] = qcins_state->x[17];
	return true;
}

bool qcins_get_thrust(uintptr_t qcins_handle, float thrust[1])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	thrust[0] = qcins_state->x[18];
	return true;
}

/******** Private methods begin here **********/

void normalize_state(float *x)
{
	float qmag = sqrtf(x[4]*x[4] + x[5]*x[5] + x[6]*x[6] + x[7]*x[7]);
	x[4] /= qmag;
	x[5] /= qmag;
	x[6] /= qmag;
	x[7] /= qmag;
}

void update_state(const float * restrict x, const float *restrict P, float *restrict xnew, float *restrict Pnew)
{
	for (int i = 0; i < NUMX; i++)
		xnew[i] = x[i];
	for (int i = 0; i < NUMP; i++)
		Pnew[i] = P[i];
}

// Restrict a change in state to a change in yaw
void update_state_z(const float * restrict xold, float *restrict xnew)
{
	float delta_yaw_direction[4] = {-xold[4+3], -xold[4+2], xold[4+1], xold[4+0]};
	float delta_q[4] = {xnew[4]-xold[4], xnew[5]-xold[5], xnew[6]-xold[6], xnew[7]-xold[7]};
	float len = delta_yaw_direction[0]*delta_q[0] + delta_yaw_direction[1]*delta_q[1] +
				delta_yaw_direction[2]*delta_q[2] + delta_yaw_direction[3]*delta_q[3];
	xnew[4+0] = xold[4+0] + delta_yaw_direction[0] * len;
	xnew[4+1] = xold[4+1] + delta_yaw_direction[1] * len;
	xnew[4+2] = xold[4+2] + delta_yaw_direction[2] * len;
	xnew[4+3] = xold[4+3] + delta_yaw_direction[3] * len;
}

void state_prediction(const float * restrict x, const float *u, float Ts, const float *restrict param, float * restrict xnew)
{
   const float x0 = 2*param[0]*x[14]*x[18];
   const float x1 = x[4]*x[6];
   const float x2 = x[5]*x[7];
   const float x3 = 2*x[4]*x[7];
   const float x4 = x3 - 2*x[5]*x[6];
   const float x5 = 2*x[3];
   const float x6 = x[4]*x[5];
   const float x7 = x[6]*x[7];
   const float x8 = x6 + x7;
   const float x9 = (x[4]*x[4]);
   const float x10 = (x[5]*x[5]);
   const float x11 = -x10 + x9;
   const float x12 = (x[6]*x[6]);
   const float x13 = (x[7]*x[7]);
   const float x14 = -x13;
   const float x15 = x11 + x12 + x14;
   const float x16 = x15*x[2] - x4*x[1] + x5*x8;
   const float x17 = x16*param[6];
   const float x18 = -x12;
   const float x19 = x10 + x14 + x18 + x9;
   const float x20 = x3 + 2*x[5]*x[6];
   const float x21 = x1 - x2;
   const float x22 = x19*x[1] + x20*x[2] - x21*x5;
   const float x23 = x22*param[6];
   const float x24 = 2*param[6];
   const float x25 = (1.0F/2.0F)*x[5];
   const float x26 = (1.0F/2.0F)*x[6];
   const float x27 = (1.0F/2.0F)*x[7];
   const float x28 = (1.0F/2.0F)*x[4];
   const float x29 = Ts/param[5];
   xnew[0] = Ts*x[3] + x[0];
   xnew[1] = Ts*(-x0*(x1 + x2) + x17*x4 - x19*x23) + x[1];
   xnew[2] = Ts*(x0*(x6 - x7) - x15*x17 - x20*x23) + x[2];
   xnew[3] = Ts*(-x16*x24*x8 + x21*x22*x24 - (x11 + x13 + x18)*param[0]*x[14]*x[18] + param[0]) + x[3];
   xnew[4] = Ts*(-x25*x[8] - x26*x[9] - x27*x[10]) + x[4];
   xnew[5] = Ts*(x26*x[10] - x27*x[9] + x28*x[8]) + x[5];
   xnew[6] = Ts*(-x25*x[10] + x27*x[8] + x28*x[9]) + x[6];
   xnew[7] = Ts*(x25*x[9] - x26*x[8] + x28*x[10]) + x[7];
   xnew[8] = Ts*param[1]*x[11] + x[8];
   xnew[9] = Ts*param[2]*x[12] + x[9];
   xnew[10] = Ts*(-(-u[2] + x[13] + x[17])*param[4] + param[3]*x[13]) + x[10];
   xnew[11] = x29*(u[0] - x[11] - x[15]) + x[11];
   xnew[12] = x29*(u[1] - x[12] - x[16]) + x[12];
   xnew[13] = x29*(u[2] - x[13] - x[17]) + x[13];
   xnew[14] = x29*(u[3] - x[14]) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = x[18];
}

void linearize_FH(const float * restrict x, const float Ts, const float * restrict param, float * restrict F, float * restrict H)
{
   const float x0 = 2*x[4];
   const float x1 = x0*x[7];
   const float x2 = 2*x[5];
   const float x3 = x2*x[6];
   const float x4 = -x1 + x3;
   const float x5 = (x[4]*x[4]);
   const float x6 = (x[5]*x[5]);
   const float x7 = (x[6]*x[6]);
   const float x8 = -x7;
   const float x9 = (x[7]*x[7]);
   const float x10 = -x9;
   const float x11 = x10 + x5 + x6 + x8;
   const float x12 = x4*param[6];
   const float x13 = x5 - x6;
   const float x14 = x10 + x13 + x7;
   const float x15 = x1 + x3;
   const float x16 = x15*param[6];
   const float x17 = Ts*(-x11*x16 - x12*x14);
   const float x18 = x0*x[5];
   const float x19 = 2*x[6];
   const float x20 = x19*x[7];
   const float x21 = x18 + x20;
   const float x22 = x21*param[6];
   const float x23 = x0*x[6];
   const float x24 = x2*x[7];
   const float x25 = -x23 + x24;
   const float x26 = x25*param[6];
   const float x27 = Ts*(-x11*x26 - x22*x4);
   const float x28 = 2*x[1];
   const float x29 = 2*x[2];
   const float x30 = 2*x[3];
   const float x31 = -x28*x[7] + x29*x[4] + x30*x[5];
   const float x32 = x28*x[4];
   const float x33 = x29*x[7];
   const float x34 = x30*x[6];
   const float x35 = x32 + x33 - x34;
   const float x36 = x11*param[6];
   const float x37 = param[0]*x[14]*x[18];
   const float x38 = x19*x37;
   const float x39 = 2*x[7];
   const float x40 = (x14*x[2] + x21*x[3] + x4*x[1])*param[6];
   const float x41 = x39*x40;
   const float x42 = (x11*x[1] + x15*x[2] + x25*x[3])*param[6];
   const float x43 = x0*x42;
   const float x44 = -x38 + x41 - x43;
   const float x45 = x28*x[6];
   const float x46 = x29*x[5];
   const float x47 = x30*x[4];
   const float x48 = x45 - x46 + x47;
   const float x49 = x28*x[5] + x29*x[6] + x30*x[7];
   const float x50 = -x19*x40 - x2*x42 - x37*x39;
   const float x51 = -x45 + x46 - x47;
   const float x52 = x0*x37;
   const float x53 = x2*x40;
   const float x54 = x19*x42;
   const float x55 = -x52 - x53 + x54;
   const float x56 = x2*x37;
   const float x57 = -x32 - x33 + x34;
   const float x58 = x0*x40;
   const float x59 = x39*x42;
   const float x60 = param[0]*x[18];
   const float x61 = Ts*(x23 + x24);
   const float x62 = param[0]*x[14];
   const float x63 = Ts*(-x14*x22 - x15*x26);
   const float x64 = x14*param[6];
   const float x65 = x56 - x58 - x59;
   const float x66 = Ts*(-x18 + x20);
   const float x67 = Ts*(x13 + x8 + x9);
   const float x68 = (1.0F/2.0F)*Ts;
   const float x69 = x68*x[8];
   const float x70 = -x69;
   const float x71 = x68*x[9];
   const float x72 = -x71;
   const float x73 = x68*x[10];
   const float x74 = -x73;
   const float x75 = x68*x[5];
   const float x76 = -x75;
   const float x77 = x68*x[6];
   const float x78 = -x77;
   const float x79 = x68*x[7];
   const float x80 = -x79;
   const float x81 = x68*x[4];
   const float x82 = -Ts/param[5];
   const float x83 = -x49*param[6];
   const float x84 = -x31*param[6];
   const float x85 = 2*param[7];
   const float x86 = x85*x[4];
   const float x87 = 2*param[8];
   const float x88 = x87*x[7];
   const float x89 = 2*param[9];
   const float x90 = x89*x[6];
   const float x91 = x86 + x88 - x90;
   const float x92 = x85*x[5] + x87*x[6] + x89*x[7];
   const float x93 = x85*x[6];
   const float x94 = x87*x[5];
   const float x95 = x89*x[4];
   const float x96 = x85*x[7];
   const float x97 = x87*x[4];
   const float x98 = x89*x[5];
   const float x99 = -x96 + x97 + x98;
   const float x100 = x93 - x94 + x95;
   F[0] = Ts;
   F[1] = Ts*(-(x11*x11)*param[6] - (x4*x4)*param[6]);
   F[2] = x17;
   F[3] = x27;
   F[4] = Ts*(-x12*x31 - x35*x36 + x44);
   F[5] = Ts*(-x12*x48 - x36*x49 + x50);
   F[6] = Ts*(-x12*x49 - x36*x51 + x55);
   F[7] = Ts*(-x12*x57 - x31*x36 - x56 + x58 + x59);
   F[8] = -x60*x61;
   F[9] = -x61*x62;
   F[10] = x17;
   F[11] = Ts*(-(x14*x14)*param[6] - (x15*x15)*param[6]);
   F[12] = x63;
   F[13] = Ts*(-x16*x35 - x31*x64 + x65);
   F[14] = Ts*(-x16*x49 - x48*x64 + x52 + x53 - x54);
   F[15] = Ts*(-x16*x51 - x49*x64 + x50);
   F[16] = Ts*(-x16*x31 + x44 - x57*x64);
   F[17] = -x60*x66;
   F[18] = -x62*x66;
   F[19] = x27;
   F[20] = x63;
   F[21] = Ts*(-(x21*x21)*param[6] - (x25*x25)*param[6]);
   F[22] = Ts*(-x22*x31 - x26*x35 + x55);
   F[23] = Ts*(-x22*x48 - x26*x49 + x65);
   F[24] = Ts*(-x22*x49 - x26*x51 + x38 - x41 + x43);
   F[25] = Ts*(-x22*x57 - x26*x31 + x50);
   F[26] = -x60*x67;
   F[27] = -x62*x67;
   F[28] = x70;
   F[29] = x72;
   F[30] = x74;
   F[31] = x76;
   F[32] = x78;
   F[33] = x80;
   F[34] = x69;
   F[35] = x73;
   F[36] = x72;
   F[37] = x81;
   F[38] = x80;
   F[39] = x77;
   F[40] = x71;
   F[41] = x74;
   F[42] = x69;
   F[43] = x79;
   F[44] = x81;
   F[45] = x76;
   F[46] = x73;
   F[47] = x71;
   F[48] = x70;
   F[49] = x78;
   F[50] = x75;
   F[51] = x81;
   F[52] = Ts*param[1];
   F[53] = Ts*param[2];
   F[54] = Ts*(param[3] - param[4]);
   F[55] = -Ts*param[4];
   F[56] = x82;
   F[57] = x82;
   F[58] = x82;
   F[59] = x82;
   F[60] = x82;
   F[61] = x82;
   F[62] = x82;
   H[0] = 1;
   H[1] = -x36;
   H[2] = -x16;
   H[3] = -x26;
   H[4] = -x35*param[6];
   H[5] = x83;
   H[6] = -x51*param[6];
   H[7] = x84;
   H[8] = -x12;
   H[9] = -x64;
   H[10] = -x22;
   H[11] = x84;
   H[12] = -x48*param[6];
   H[13] = x83;
   H[14] = -x57*param[6];
   H[15] = -x60;
   H[16] = -x62;
   H[17] = 1;
   H[18] = 1;
   H[19] = 1;
   H[20] = x91;
   H[21] = x92;
   H[22] = -x93 + x94 - x95;
   H[23] = x99;
   H[24] = x99;
   H[25] = x100;
   H[26] = x92;
   H[27] = -x86 - x88 + x90;
   H[28] = x100;
   H[29] = x96 - x97 - x98;
   H[30] = x91;
   H[31] = x92;
}

// TODO: rework F definition so that Ts is already multiplied into all the terms
// to reduce the memory requirement (currently duplicates F)
void covariance_prediction(const float *restrict P, const float * restrict F, const float * restrict Q, float * restrict Pnew)
{
   const float x0 = F[0]*P[36] + P[3];
   const float x1 = F[0]*P[26] + P[2];
   const float x2 = F[0]*P[37] + P[4];
   const float x3 = F[0]*P[38] + P[5];
   const float x4 = F[0]*P[39] + P[6];
   const float x5 = F[0]*P[40] + P[7];
   const float x6 = F[0]*P[44] + P[11];
   const float x7 = F[0]*P[45] + P[12];
   const float x8 = F[1] + 1;
   const float x9 = F[0]*P[15] + P[1];
   const float x10 = F[11] + 1;
   const float x11 = F[21] + 1;
   const float x12 = F[0]*P[41] + P[8];
   const float x13 = F[0]*P[42] + P[9];
   const float x14 = F[0]*P[43] + P[10];
   const float x15 = F[62] + 1;
   const float x16 = x8*P[14] + F[2]*P[25] + F[3]*P[26] + F[4]*P[27] + F[5]*P[28] + F[6]*P[29] + F[7]*P[30] + F[8]*P[34] + F[9]*P[35];
   const float x17 = x8*P[15] + F[2]*P[26] + F[3]*P[36] + F[4]*P[37] + F[5]*P[38] + F[6]*P[39] + F[7]*P[40] + F[8]*P[44] + F[9]*P[45];
   const float x18 = x8*P[16] + F[2]*P[27] + F[3]*P[37] + F[4]*P[46] + F[5]*P[47] + F[6]*P[48] + F[7]*P[49] + F[8]*P[53] + F[9]*P[54];
   const float x19 = x8*P[17] + F[2]*P[28] + F[3]*P[38] + F[4]*P[47] + F[5]*P[55] + F[6]*P[56] + F[7]*P[57] + F[8]*P[61] + F[9]*P[62];
   const float x20 = x8*P[18] + F[2]*P[29] + F[3]*P[39] + F[4]*P[48] + F[5]*P[56] + F[6]*P[63] + F[7]*P[64] + F[8]*P[68] + F[9]*P[69];
   const float x21 = x8*P[19] + F[2]*P[30] + F[3]*P[40] + F[4]*P[49] + F[5]*P[57] + F[6]*P[64] + F[7]*P[70] + F[8]*P[74] + F[9]*P[75];
   const float x22 = x8*P[23] + F[2]*P[34] + F[3]*P[44] + F[4]*P[53] + F[5]*P[61] + F[6]*P[68] + F[7]*P[74] + F[8]*P[91] + F[9]*P[92];
   const float x23 = x8*P[24] + F[2]*P[35] + F[3]*P[45] + F[4]*P[54] + F[5]*P[62] + F[6]*P[69] + F[7]*P[75] + F[8]*P[92] + F[9]*P[96];
   const float x24 = x8*P[13] + F[2]*P[14] + F[3]*P[15] + F[4]*P[16] + F[5]*P[17] + F[6]*P[18] + F[7]*P[19] + F[8]*P[23] + F[9]*P[24];
   const float x25 = x8*P[20] + F[2]*P[31] + F[3]*P[41] + F[4]*P[50] + F[5]*P[58] + F[6]*P[65] + F[7]*P[71];
   const float x26 = x8*P[21] + F[2]*P[32] + F[3]*P[42] + F[4]*P[51] + F[5]*P[59] + F[6]*P[66] + F[7]*P[72];
   const float x27 = x8*P[22] + F[2]*P[33] + F[3]*P[43] + F[4]*P[52] + F[5]*P[60] + F[6]*P[67] + F[7]*P[73];
   const float x28 = x10*P[14] + F[10]*P[13] + F[12]*P[15] + F[13]*P[16] + F[14]*P[17] + F[15]*P[18] + F[16]*P[19] + F[17]*P[23] + F[18]*P[24];
   const float x29 = x10*P[26] + F[10]*P[15] + F[12]*P[36] + F[13]*P[37] + F[14]*P[38] + F[15]*P[39] + F[16]*P[40] + F[17]*P[44] + F[18]*P[45];
   const float x30 = x10*P[27] + F[10]*P[16] + F[12]*P[37] + F[13]*P[46] + F[14]*P[47] + F[15]*P[48] + F[16]*P[49] + F[17]*P[53] + F[18]*P[54];
   const float x31 = x10*P[28] + F[10]*P[17] + F[12]*P[38] + F[13]*P[47] + F[14]*P[55] + F[15]*P[56] + F[16]*P[57] + F[17]*P[61] + F[18]*P[62];
   const float x32 = x10*P[29] + F[10]*P[18] + F[12]*P[39] + F[13]*P[48] + F[14]*P[56] + F[15]*P[63] + F[16]*P[64] + F[17]*P[68] + F[18]*P[69];
   const float x33 = x10*P[30] + F[10]*P[19] + F[12]*P[40] + F[13]*P[49] + F[14]*P[57] + F[15]*P[64] + F[16]*P[70] + F[17]*P[74] + F[18]*P[75];
   const float x34 = x10*P[34] + F[10]*P[23] + F[12]*P[44] + F[13]*P[53] + F[14]*P[61] + F[15]*P[68] + F[16]*P[74] + F[17]*P[91] + F[18]*P[92];
   const float x35 = x10*P[35] + F[10]*P[24] + F[12]*P[45] + F[13]*P[54] + F[14]*P[62] + F[15]*P[69] + F[16]*P[75] + F[17]*P[92] + F[18]*P[96];
   const float x36 = x10*P[25] + F[10]*P[14] + F[12]*P[26] + F[13]*P[27] + F[14]*P[28] + F[15]*P[29] + F[16]*P[30] + F[17]*P[34] + F[18]*P[35];
   const float x37 = x10*P[31] + F[10]*P[20] + F[12]*P[41] + F[13]*P[50] + F[14]*P[58] + F[15]*P[65] + F[16]*P[71];
   const float x38 = x10*P[32] + F[10]*P[21] + F[12]*P[42] + F[13]*P[51] + F[14]*P[59] + F[15]*P[66] + F[16]*P[72];
   const float x39 = x10*P[33] + F[10]*P[22] + F[12]*P[43] + F[13]*P[52] + F[14]*P[60] + F[15]*P[67] + F[16]*P[73];
   const float x40 = x11*P[37] + F[19]*P[16] + F[20]*P[27] + F[22]*P[46] + F[23]*P[47] + F[24]*P[48] + F[25]*P[49] + F[26]*P[53] + F[27]*P[54];
   const float x41 = x11*P[38] + F[19]*P[17] + F[20]*P[28] + F[22]*P[47] + F[23]*P[55] + F[24]*P[56] + F[25]*P[57] + F[26]*P[61] + F[27]*P[62];
   const float x42 = x11*P[39] + F[19]*P[18] + F[20]*P[29] + F[22]*P[48] + F[23]*P[56] + F[24]*P[63] + F[25]*P[64] + F[26]*P[68] + F[27]*P[69];
   const float x43 = x11*P[40] + F[19]*P[19] + F[20]*P[30] + F[22]*P[49] + F[23]*P[57] + F[24]*P[64] + F[25]*P[70] + F[26]*P[74] + F[27]*P[75];
   const float x44 = x11*P[44] + F[19]*P[23] + F[20]*P[34] + F[22]*P[53] + F[23]*P[61] + F[24]*P[68] + F[25]*P[74] + F[26]*P[91] + F[27]*P[92];
   const float x45 = x11*P[45] + F[19]*P[24] + F[20]*P[35] + F[22]*P[54] + F[23]*P[62] + F[24]*P[69] + F[25]*P[75] + F[26]*P[92] + F[27]*P[96];
   const float x46 = x11*P[41] + F[19]*P[20] + F[20]*P[31] + F[22]*P[50] + F[23]*P[58] + F[24]*P[65] + F[25]*P[71];
   const float x47 = x11*P[42] + F[19]*P[21] + F[20]*P[32] + F[22]*P[51] + F[23]*P[59] + F[24]*P[66] + F[25]*P[72];
   const float x48 = x11*P[43] + F[19]*P[22] + F[20]*P[33] + F[22]*P[52] + F[23]*P[60] + F[24]*P[67] + F[25]*P[73];
   const float x49 = F[28]*P[47] + F[29]*P[48] + F[30]*P[49] + F[31]*P[50] + F[32]*P[51] + F[33]*P[52] + P[46];
   const float x50 = F[28]*P[58] + F[29]*P[65] + F[30]*P[71] + F[31]*P[76] + P[50];
   const float x51 = F[28]*P[59] + F[29]*P[66] + F[30]*P[72] + F[32]*P[80] + P[51];
   const float x52 = F[28]*P[60] + F[29]*P[67] + F[30]*P[73] + F[33]*P[82] + P[52];
   const float x53 = F[28]*P[55] + F[29]*P[56] + F[30]*P[57] + F[31]*P[58] + F[32]*P[59] + F[33]*P[60] + P[47];
   const float x54 = F[28]*P[56] + F[29]*P[63] + F[30]*P[64] + F[31]*P[65] + F[32]*P[66] + F[33]*P[67] + P[48];
   const float x55 = F[28]*P[57] + F[29]*P[64] + F[30]*P[70] + F[31]*P[71] + F[32]*P[72] + F[33]*P[73] + P[49];
   const float x56 = F[52]*P[77];
   const float x57 = F[53]*P[78];
   const float x58 = F[54]*P[83];
   const float x59 = F[55]*P[84];
   const float x60 = F[34]*P[47] + F[35]*P[56] + F[36]*P[57] + F[37]*P[58] + F[38]*P[59] + F[39]*P[60] + P[55];
   const float x61 = F[34]*P[50] + F[35]*P[65] + F[36]*P[71] + F[37]*P[76] + P[58];
   const float x62 = F[34]*P[51] + F[35]*P[66] + F[36]*P[72] + F[38]*P[80] + P[59];
   const float x63 = F[34]*P[52] + F[35]*P[67] + F[36]*P[73] + F[39]*P[82] + P[60];
   const float x64 = F[34]*P[46] + F[35]*P[48] + F[36]*P[49] + F[37]*P[50] + F[38]*P[51] + F[39]*P[52] + P[47];
   const float x65 = F[34]*P[48] + F[35]*P[63] + F[36]*P[64] + F[37]*P[65] + F[38]*P[66] + F[39]*P[67] + P[56];
   const float x66 = F[34]*P[49] + F[35]*P[64] + F[36]*P[70] + F[37]*P[71] + F[38]*P[72] + F[39]*P[73] + P[57];
   const float x67 = F[40]*P[48] + F[41]*P[56] + F[42]*P[64] + F[43]*P[65] + F[44]*P[66] + F[45]*P[67] + P[63];
   const float x68 = F[40]*P[50] + F[41]*P[58] + F[42]*P[71] + F[43]*P[76] + P[65];
   const float x69 = F[40]*P[51] + F[41]*P[59] + F[42]*P[72] + F[44]*P[80] + P[66];
   const float x70 = F[40]*P[52] + F[41]*P[60] + F[42]*P[73] + F[45]*P[82] + P[67];
   const float x71 = F[40]*P[46] + F[41]*P[47] + F[42]*P[49] + F[43]*P[50] + F[44]*P[51] + F[45]*P[52] + P[48];
   const float x72 = F[40]*P[47] + F[41]*P[55] + F[42]*P[57] + F[43]*P[58] + F[44]*P[59] + F[45]*P[60] + P[56];
   const float x73 = F[40]*P[49] + F[41]*P[57] + F[42]*P[70] + F[43]*P[71] + F[44]*P[72] + F[45]*P[73] + P[64];
   const float x74 = F[46]*P[50] + F[47]*P[58] + F[48]*P[65] + F[49]*P[76] + P[71];
   const float x75 = F[46]*P[51] + F[47]*P[59] + F[48]*P[66] + F[50]*P[80] + P[72];
   const float x76 = F[46]*P[52] + F[47]*P[60] + F[48]*P[67] + F[51]*P[82] + P[73];
   const float x77 = F[52]*P[85] + P[77];
   const float x78 = F[52]*P[86] + P[79];
   const float x79 = F[56] + 1;
   const float x80 = F[58] + 1;
   const float x81 = F[54]*P[89] + F[55]*P[90] + P[83];
   const float x82 = F[54]*P[90] + F[55]*P[95] + P[84];
   const float x83 = F[60] + 1;
   const float x84 = x79*P[86] + F[57]*P[93];
   const float x85 = x80*P[88] + F[59]*P[94];
   const float x86 = x83*P[90] + F[61]*P[95];
   Pnew[0] = x0*F[0] + F[0]*P[3] + P[0] + Q[0];
   Pnew[1] = x0*F[3] + x1*F[2] + x2*F[4] + x3*F[5] + x4*F[6] + x5*F[7] + x6*F[8] + x7*F[9] + x8*x9;
   Pnew[2] = x0*F[12] + x1*x10 + x2*F[13] + x3*F[14] + x4*F[15] + x5*F[16] + x6*F[17] + x7*F[18] + x9*F[10];
   Pnew[3] = x0*x11 + x1*F[20] + x2*F[22] + x3*F[23] + x4*F[24] + x5*F[25] + x6*F[26] + x7*F[27] + x9*F[19];
   Pnew[4] = x12*F[31] + x13*F[32] + x14*F[33] + x2 + x3*F[28] + x4*F[29] + x5*F[30];
   Pnew[5] = x12*F[37] + x13*F[38] + x14*F[39] + x2*F[34] + x3 + x4*F[35] + x5*F[36];
   Pnew[6] = x12*F[43] + x13*F[44] + x14*F[45] + x2*F[40] + x3*F[41] + x4 + x5*F[42];
   Pnew[7] = x12*F[49] + x13*F[50] + x14*F[51] + x2*F[46] + x3*F[47] + x4*F[48] + x5;
   Pnew[8] = x12;
   Pnew[9] = x13;
   Pnew[10] = x14;
   Pnew[11] = x15*x6;
   Pnew[12] = x7;
   Pnew[13] = x16*F[2] + x17*F[3] + x18*F[4] + x19*F[5] + x20*F[6] + x21*F[7] + x22*F[8] + x23*F[9] + x24*x8 + Q[1];
   Pnew[14] = x10*x16 + x17*F[12] + x18*F[13] + x19*F[14] + x20*F[15] + x21*F[16] + x22*F[17] + x23*F[18] + x24*F[10];
   Pnew[15] = x11*x17 + x16*F[20] + x18*F[22] + x19*F[23] + x20*F[24] + x21*F[25] + x22*F[26] + x23*F[27] + x24*F[19];
   Pnew[16] = x18 + x19*F[28] + x20*F[29] + x21*F[30] + x25*F[31] + x26*F[32] + x27*F[33];
   Pnew[17] = x18*F[34] + x19 + x20*F[35] + x21*F[36] + x25*F[37] + x26*F[38] + x27*F[39];
   Pnew[18] = x18*F[40] + x19*F[41] + x20 + x21*F[42] + x25*F[43] + x26*F[44] + x27*F[45];
   Pnew[19] = x18*F[46] + x19*F[47] + x20*F[48] + x21 + x25*F[49] + x26*F[50] + x27*F[51];
   Pnew[20] = x25;
   Pnew[21] = x26;
   Pnew[22] = x27;
   Pnew[23] = x15*x22;
   Pnew[24] = x23;
   Pnew[25] = x10*x36 + x28*F[10] + x29*F[12] + x30*F[13] + x31*F[14] + x32*F[15] + x33*F[16] + x34*F[17] + x35*F[18] + Q[2];
   Pnew[26] = x11*x29 + x28*F[19] + x30*F[22] + x31*F[23] + x32*F[24] + x33*F[25] + x34*F[26] + x35*F[27] + x36*F[20];
   Pnew[27] = x30 + x31*F[28] + x32*F[29] + x33*F[30] + x37*F[31] + x38*F[32] + x39*F[33];
   Pnew[28] = x30*F[34] + x31 + x32*F[35] + x33*F[36] + x37*F[37] + x38*F[38] + x39*F[39];
   Pnew[29] = x30*F[40] + x31*F[41] + x32 + x33*F[42] + x37*F[43] + x38*F[44] + x39*F[45];
   Pnew[30] = x30*F[46] + x31*F[47] + x32*F[48] + x33 + x37*F[49] + x38*F[50] + x39*F[51];
   Pnew[31] = x37;
   Pnew[32] = x38;
   Pnew[33] = x39;
   Pnew[34] = x15*x34;
   Pnew[35] = x35;
   Pnew[36] = x11*(x11*P[36] + F[19]*P[15] + F[20]*P[26] + F[22]*P[37] + F[23]*P[38] + F[24]*P[39] + F[25]*P[40] + F[26]*P[44] + F[27]*P[45]) + x40*F[22] + x41*F[23] + x42*F[24] + x43*F[25] + x44*F[26] + x45*F[27] + (x11*P[15] + F[19]*P[13] + F[20]*P[14] + F[22]*P[16] + F[23]*P[17] + F[24]*P[18] + F[25]*P[19] + F[26]*P[23] + F[27]*P[24])*F[19] + (x11*P[26] + F[19]*P[14] + F[20]*P[25] + F[22]*P[27] + F[23]*P[28] + F[24]*P[29] + F[25]*P[30] + F[26]*P[34] + F[27]*P[35])*F[20] + Q[3];
   Pnew[37] = x40 + x41*F[28] + x42*F[29] + x43*F[30] + x46*F[31] + x47*F[32] + x48*F[33];
   Pnew[38] = x40*F[34] + x41 + x42*F[35] + x43*F[36] + x46*F[37] + x47*F[38] + x48*F[39];
   Pnew[39] = x40*F[40] + x41*F[41] + x42 + x43*F[42] + x46*F[43] + x47*F[44] + x48*F[45];
   Pnew[40] = x40*F[46] + x41*F[47] + x42*F[48] + x43 + x46*F[49] + x47*F[50] + x48*F[51];
   Pnew[41] = x46;
   Pnew[42] = x47;
   Pnew[43] = x48;
   Pnew[44] = x15*x44;
   Pnew[45] = x45;
   Pnew[46] = x49 + x50*F[31] + x51*F[32] + x52*F[33] + x53*F[28] + x54*F[29] + x55*F[30] + Q[4];
   Pnew[47] = x49*F[34] + x50*F[37] + x51*F[38] + x52*F[39] + x53 + x54*F[35] + x55*F[36];
   Pnew[48] = x49*F[40] + x50*F[43] + x51*F[44] + x52*F[45] + x53*F[41] + x54 + x55*F[42];
   Pnew[49] = x49*F[46] + x50*F[49] + x51*F[50] + x52*F[51] + x53*F[47] + x54*F[48] + x55;
   Pnew[50] = x50 + x56*F[31];
   Pnew[51] = x51 + x57*F[31];
   Pnew[52] = x52 + x58*F[33] + x59*F[33];
   Pnew[53] = x15*(F[28]*P[61] + F[29]*P[68] + F[30]*P[74] + P[53]);
   Pnew[54] = F[28]*P[62] + F[29]*P[69] + F[30]*P[75] + P[54];
   Pnew[55] = x60 + x61*F[37] + x62*F[38] + x63*F[39] + x64*F[34] + x65*F[35] + x66*F[36] + Q[5];
   Pnew[56] = x60*F[41] + x61*F[43] + x62*F[44] + x63*F[45] + x64*F[40] + x65 + x66*F[42];
   Pnew[57] = x60*F[47] + x61*F[49] + x62*F[50] + x63*F[51] + x64*F[46] + x65*F[48] + x66;
   Pnew[58] = x56*F[37] + x61;
   Pnew[59] = x57*F[37] + x62;
   Pnew[60] = x58*F[39] + x59*F[39] + x63;
   Pnew[61] = x15*(F[34]*P[53] + F[35]*P[68] + F[36]*P[74] + P[61]);
   Pnew[62] = F[34]*P[54] + F[35]*P[69] + F[36]*P[75] + P[62];
   Pnew[63] = x67 + x68*F[43] + x69*F[44] + x70*F[45] + x71*F[40] + x72*F[41] + x73*F[42] + Q[6];
   Pnew[64] = x67*F[48] + x68*F[49] + x69*F[50] + x70*F[51] + x71*F[46] + x72*F[47] + x73;
   Pnew[65] = x56*F[43] + x68;
   Pnew[66] = x57*F[43] + x69;
   Pnew[67] = x58*F[45] + x59*F[45] + x70;
   Pnew[68] = x15*(F[40]*P[53] + F[41]*P[61] + F[42]*P[74] + P[68]);
   Pnew[69] = F[40]*P[54] + F[41]*P[62] + F[42]*P[75] + P[69];
   Pnew[70] = x74*F[49] + x75*F[50] + x76*F[51] + (F[46]*P[46] + F[47]*P[47] + F[48]*P[48] + F[49]*P[50] + F[50]*P[51] + F[51]*P[52] + P[49])*F[46] + (F[46]*P[47] + F[47]*P[55] + F[48]*P[56] + F[49]*P[58] + F[50]*P[59] + F[51]*P[60] + P[57])*F[47] + (F[46]*P[48] + F[47]*P[56] + F[48]*P[63] + F[49]*P[65] + F[50]*P[66] + F[51]*P[67] + P[64])*F[48] + F[46]*P[49] + F[47]*P[57] + F[48]*P[64] + F[49]*P[71] + F[50]*P[72] + F[51]*P[73] + P[70] + Q[7];
   Pnew[71] = x56*F[49] + x74;
   Pnew[72] = x57*F[49] + x75;
   Pnew[73] = x58*F[51] + x59*F[51] + x76;
   Pnew[74] = x15*(F[46]*P[53] + F[47]*P[61] + F[48]*P[68] + P[74]);
   Pnew[75] = F[46]*P[54] + F[47]*P[62] + F[48]*P[69] + P[75];
   Pnew[76] = x56 + x77*F[52] + P[76] + Q[8];
   Pnew[77] = x77*x79 + x78*F[57];
   Pnew[78] = x80*P[78];
   Pnew[79] = x78;
   Pnew[80] = (F[53]*F[53])*P[87] + P[80] + Q[9];
   Pnew[81] = F[53]*P[88] + P[81];
   Pnew[82] = x58 + x59 + x81*F[54] + x82*F[55] + P[82] + Q[10];
   Pnew[83] = x81*x83 + x82*F[61];
   Pnew[84] = x82;
   Pnew[85] = x79*(x79*P[85] + F[57]*P[86]) + x84*F[57] + Q[11];
   Pnew[86] = x84;
   Pnew[87] = x80*(x80*P[87] + F[59]*P[88]) + x85*F[59] + Q[12];
   Pnew[88] = x85;
   Pnew[89] = x83*(x83*P[89] + F[61]*P[90]) + x86*F[61] + Q[13];
   Pnew[90] = x86;
   Pnew[91] = (x15*x15)*P[91] + Q[14];
   Pnew[92] = x15*P[92];
   Pnew[93] = P[93] + Q[15];
   Pnew[94] = P[94] + Q[16];
   Pnew[95] = P[95] + Q[17];
   Pnew[96] = P[96] + Q[18];
}

void baro_correction(const float *restrict x, const float *restrict P, float baro, const float *restrict H,
	                      const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = (H[0]*H[0]);
   const float x1 = x0*P[0];
   const float x2 = 1.0F/(x1 + R[0]);
   const float x3 = x2*(-baro + x[0])*H[0];
   const float x4 = -x1*x2 + 1;
   const float x5 = x0*x2;
   const float x6 = x0*x2*P[1];
   const float x7 = x0*x2*P[2];
   const float x8 = x0*x2*P[3];
   const float x9 = x0*x2*P[4];
   const float x10 = x0*x2*P[5];
   const float x11 = x0*x2*P[6];
   const float x12 = x0*x2*P[7];
   xnew[0] = -x3*P[0] + x[0];
   xnew[1] = -x3*P[1] + x[1];
   xnew[2] = -x3*P[2] + x[2];
   xnew[3] = -x3*P[3] + x[3];
   xnew[4] = -x3*P[4] + x[4];
   xnew[5] = -x3*P[5] + x[5];
   xnew[6] = -x3*P[6] + x[6];
   xnew[7] = -x3*P[7] + x[7];
   xnew[8] = -x3*P[8] + x[8];
   xnew[9] = -x3*P[9] + x[9];
   xnew[10] = -x3*P[10] + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = -x3*P[11] + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = -x3*P[12] + x[18];
   Pnew[0] = x4*P[0];
   Pnew[1] = x4*P[1];
   Pnew[2] = x4*P[2];
   Pnew[3] = x4*P[3];
   Pnew[4] = x4*P[4];
   Pnew[5] = x4*P[5];
   Pnew[6] = x4*P[6];
   Pnew[7] = x4*P[7];
   Pnew[8] = x4*P[8];
   Pnew[9] = x4*P[9];
   Pnew[10] = x4*P[10];
   Pnew[11] = x4*P[11];
   Pnew[12] = x4*P[12];
   Pnew[13] = -x5*(P[1]*P[1]) + P[13];
   Pnew[14] = -x6*P[2] + P[14];
   Pnew[15] = -x6*P[3] + P[15];
   Pnew[16] = -x6*P[4] + P[16];
   Pnew[17] = -x6*P[5] + P[17];
   Pnew[18] = -x6*P[6] + P[18];
   Pnew[19] = -x6*P[7] + P[19];
   Pnew[20] = -x6*P[8] + P[20];
   Pnew[21] = -x6*P[9] + P[21];
   Pnew[22] = -x6*P[10] + P[22];
   Pnew[23] = -x6*P[11] + P[23];
   Pnew[24] = -x6*P[12] + P[24];
   Pnew[25] = -x5*(P[2]*P[2]) + P[25];
   Pnew[26] = -x7*P[3] + P[26];
   Pnew[27] = -x7*P[4] + P[27];
   Pnew[28] = -x7*P[5] + P[28];
   Pnew[29] = -x7*P[6] + P[29];
   Pnew[30] = -x7*P[7] + P[30];
   Pnew[31] = -x7*P[8] + P[31];
   Pnew[32] = -x7*P[9] + P[32];
   Pnew[33] = -x7*P[10] + P[33];
   Pnew[34] = -x7*P[11] + P[34];
   Pnew[35] = -x7*P[12] + P[35];
   Pnew[36] = -x5*(P[3]*P[3]) + P[36];
   Pnew[37] = -x8*P[4] + P[37];
   Pnew[38] = -x8*P[5] + P[38];
   Pnew[39] = -x8*P[6] + P[39];
   Pnew[40] = -x8*P[7] + P[40];
   Pnew[41] = -x8*P[8] + P[41];
   Pnew[42] = -x8*P[9] + P[42];
   Pnew[43] = -x8*P[10] + P[43];
   Pnew[44] = -x8*P[11] + P[44];
   Pnew[45] = -x8*P[12] + P[45];
   Pnew[46] = -x5*(P[4]*P[4]) + P[46];
   Pnew[47] = -x9*P[5] + P[47];
   Pnew[48] = -x9*P[6] + P[48];
   Pnew[49] = -x9*P[7] + P[49];
   Pnew[50] = -x9*P[8] + P[50];
   Pnew[51] = -x9*P[9] + P[51];
   Pnew[52] = -x9*P[10] + P[52];
   Pnew[53] = -x9*P[11] + P[53];
   Pnew[54] = -x9*P[12] + P[54];
   Pnew[55] = -x5*(P[5]*P[5]) + P[55];
   Pnew[56] = -x10*P[6] + P[56];
   Pnew[57] = -x10*P[7] + P[57];
   Pnew[58] = -x10*P[8] + P[58];
   Pnew[59] = -x10*P[9] + P[59];
   Pnew[60] = -x10*P[10] + P[60];
   Pnew[61] = -x10*P[11] + P[61];
   Pnew[62] = -x10*P[12] + P[62];
   Pnew[63] = -x5*(P[6]*P[6]) + P[63];
   Pnew[64] = -x11*P[7] + P[64];
   Pnew[65] = -x11*P[8] + P[65];
   Pnew[66] = -x11*P[9] + P[66];
   Pnew[67] = -x11*P[10] + P[67];
   Pnew[68] = -x11*P[11] + P[68];
   Pnew[69] = -x11*P[12] + P[69];
   Pnew[70] = -x5*(P[7]*P[7]) + P[70];
   Pnew[71] = -x12*P[8] + P[71];
   Pnew[72] = -x12*P[9] + P[72];
   Pnew[73] = -x12*P[10] + P[73];
   Pnew[74] = -x12*P[11] + P[74];
   Pnew[75] = -x12*P[12] + P[75];
   Pnew[76] = -x5*(P[8]*P[8]) + P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = -x5*(P[9]*P[9]) + P[80];
   Pnew[81] = P[81];
   Pnew[82] = -x5*(P[10]*P[10]) + P[82];
   Pnew[83] = P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
   Pnew[87] = P[87];
   Pnew[88] = P[88];
   Pnew[89] = P[89];
   Pnew[90] = P[90];
   Pnew[91] = -x5*(P[11]*P[11]) + P[91];
   Pnew[92] = -x5*P[11]*P[12] + P[92];
   Pnew[93] = P[93];
   Pnew[94] = P[94];
   Pnew[95] = P[95];
   Pnew[96] = -x5*(P[12]*P[12]) + P[96];
}

void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, const float * restrict param,
        float *restrict xnew, float *restrict Pnew)
{
   const float x0 = 2*param[7];
   const float x1 = x[4]*x[6];
   const float x2 = x[5]*x[7];
   const float x3 = 2*param[8];
   const float x4 = x[4]*x[5];
   const float x5 = x[6]*x[7];
   const float x6 = (x[4]*x[4]);
   const float x7 = (x[5]*x[5]);
   const float x8 = x6 - x7;
   const float x9 = (x[7]*x[7]);
   const float x10 = (x[6]*x[6]);
   const float x11 = -x10;
   const float x12 = x0*(x1 + x2) - x3*(x4 - x5) + (x11 + x8 + x9)*param[9] - mag[2];
   const float x13 = H[20]*P[4] + H[21]*P[5] + H[22]*P[6] + H[23]*P[7];
   const float x14 = H[20]*P[46] + H[21]*P[47] + H[22]*P[48] + H[23]*P[49];
   const float x15 = H[20]*P[47] + H[21]*P[55] + H[22]*P[56] + H[23]*P[57];
   const float x16 = H[20]*P[48] + H[21]*P[56] + H[22]*P[63] + H[23]*P[64];
   const float x17 = H[20]*P[49] + H[21]*P[57] + H[22]*P[64] + H[23]*P[70];
   const float x18 = x14*H[24] + x15*H[25] + x16*H[26] + x17*H[27];
   const float x19 = H[24]*P[46] + H[25]*P[47] + H[26]*P[48] + H[27]*P[49];
   const float x20 = H[24]*P[47] + H[25]*P[55] + H[26]*P[56] + H[27]*P[57];
   const float x21 = H[24]*P[48] + H[25]*P[56] + H[26]*P[63] + H[27]*P[64];
   const float x22 = H[24]*P[49] + H[25]*P[57] + H[26]*P[64] + H[27]*P[70];
   const float x23 = x19*H[28] + x20*H[29] + x21*H[30] + x22*H[31];
   const float x24 = x18*x23;
   const float x25 = x14*H[28] + x15*H[29] + x16*H[30] + x17*H[31];
   const float x26 = x19*H[24] + x20*H[25] + x21*H[26] + x22*H[27] + R[8];
   const float x27 = H[28]*P[46] + H[29]*P[47] + H[30]*P[48] + H[31]*P[49];
   const float x28 = H[28]*P[47] + H[29]*P[55] + H[30]*P[56] + H[31]*P[57];
   const float x29 = H[28]*P[48] + H[29]*P[56] + H[30]*P[63] + H[31]*P[64];
   const float x30 = H[28]*P[49] + H[29]*P[57] + H[30]*P[64] + H[31]*P[70];
   const float x31 = x27*H[20] + x28*H[21] + x29*H[22] + x30*H[23];
   const float x32 = x27*H[24] + x28*H[25] + x29*H[26] + x30*H[27];
   const float x33 = x19*H[20] + x20*H[21] + x21*H[22] + x22*H[23];
   const float x34 = x25*x33;
   const float x35 = x27*H[28] + x28*H[29] + x29*H[30] + x30*H[31] + R[9];
   const float x36 = x18*x33;
   const float x37 = x25*x31;
   const float x38 = x14*H[20] + x15*H[21] + x16*H[22] + x17*H[23] + R[7];
   const float x39 = x23*x32;
   const float x40 = x26*x38;
   const float x41 = 1.0F/(x24*x31 - x26*x37 + x32*x34 - x35*x36 + x35*x40 - x38*x39);
   const float x42 = x41*(x24 - x25*x26);
   const float x43 = x13*x42;
   const float x44 = H[24]*P[4] + H[25]*P[5] + H[26]*P[6] + H[27]*P[7];
   const float x45 = x41*(-x23*x38 + x34);
   const float x46 = x44*x45;
   const float x47 = H[28]*P[4] + H[29]*P[5] + H[30]*P[6] + H[31]*P[7];
   const float x48 = x41*(-x36 + x40);
   const float x49 = x47*x48;
   const float x50 = 2*param[9];
   const float x51 = x[4]*x[7];
   const float x52 = x[5]*x[6];
   const float x53 = -x9;
   const float x54 = -x0*(x51 - x52) + x50*(x4 + x5) + (x10 + x53 + x8)*param[8] - mag[1];
   const float x55 = x41*(x18*x31 - x32*x38);
   const float x56 = x47*x55;
   const float x57 = x41*(-x18*x35 + x25*x32);
   const float x58 = x13*x57;
   const float x59 = x41*(x35*x38 - x37);
   const float x60 = x44*x59;
   const float x61 = x3*(x51 + x52) - x50*(x1 - x2) + (x11 + x53 + x6 + x7)*param[7] - mag[0];
   const float x62 = x41*(-x26*x31 + x32*x33);
   const float x63 = x47*x62;
   const float x64 = x41*(x23*x31 - x33*x35);
   const float x65 = x44*x64;
   const float x66 = x41*(x26*x35 - x39);
   const float x67 = x13*x66;
   const float x68 = H[20]*P[16] + H[21]*P[17] + H[22]*P[18] + H[23]*P[19];
   const float x69 = x42*x68;
   const float x70 = H[24]*P[16] + H[25]*P[17] + H[26]*P[18] + H[27]*P[19];
   const float x71 = x45*x70;
   const float x72 = H[28]*P[16] + H[29]*P[17] + H[30]*P[18] + H[31]*P[19];
   const float x73 = x48*x72;
   const float x74 = x55*x72;
   const float x75 = x57*x68;
   const float x76 = x59*x70;
   const float x77 = x62*x72;
   const float x78 = x64*x70;
   const float x79 = x66*x68;
   const float x80 = H[20]*P[27] + H[21]*P[28] + H[22]*P[29] + H[23]*P[30];
   const float x81 = x42*x80;
   const float x82 = H[24]*P[27] + H[25]*P[28] + H[26]*P[29] + H[27]*P[30];
   const float x83 = x45*x82;
   const float x84 = H[28]*P[27] + H[29]*P[28] + H[30]*P[29] + H[31]*P[30];
   const float x85 = x48*x84;
   const float x86 = x55*x84;
   const float x87 = x57*x80;
   const float x88 = x59*x82;
   const float x89 = x62*x84;
   const float x90 = x64*x82;
   const float x91 = x66*x80;
   const float x92 = H[20]*P[37] + H[21]*P[38] + H[22]*P[39] + H[23]*P[40];
   const float x93 = x42*x92;
   const float x94 = H[24]*P[37] + H[25]*P[38] + H[26]*P[39] + H[27]*P[40];
   const float x95 = x45*x94;
   const float x96 = H[28]*P[37] + H[29]*P[38] + H[30]*P[39] + H[31]*P[40];
   const float x97 = x48*x96;
   const float x98 = x55*x96;
   const float x99 = x57*x92;
   const float x100 = x59*x94;
   const float x101 = x62*x96;
   const float x102 = x64*x94;
   const float x103 = x66*x92;
   const float x104 = x14*x42;
   const float x105 = x19*x45;
   const float x106 = x27*x48;
   const float x107 = x27*x55;
   const float x108 = x14*x57;
   const float x109 = x19*x59;
   const float x110 = x27*x62;
   const float x111 = x19*x64;
   const float x112 = x14*x66;
   const float x113 = x15*x42;
   const float x114 = x20*x45;
   const float x115 = x28*x48;
   const float x116 = x28*x55;
   const float x117 = x15*x57;
   const float x118 = x20*x59;
   const float x119 = x28*x62;
   const float x120 = x20*x64;
   const float x121 = x15*x66;
   const float x122 = x16*x42;
   const float x123 = x21*x45;
   const float x124 = x29*x48;
   const float x125 = x29*x55;
   const float x126 = x16*x57;
   const float x127 = x21*x59;
   const float x128 = x29*x62;
   const float x129 = x21*x64;
   const float x130 = x16*x66;
   const float x131 = x17*x42;
   const float x132 = x22*x45;
   const float x133 = x30*x48;
   const float x134 = x30*x55;
   const float x135 = x17*x57;
   const float x136 = x22*x59;
   const float x137 = x30*x62;
   const float x138 = x22*x64;
   const float x139 = x17*x66;
   const float x140 = H[20]*P[50] + H[21]*P[58] + H[22]*P[65] + H[23]*P[71];
   const float x141 = x140*x42;
   const float x142 = H[24]*P[50] + H[25]*P[58] + H[26]*P[65] + H[27]*P[71];
   const float x143 = x142*x45;
   const float x144 = H[28]*P[50] + H[29]*P[58] + H[30]*P[65] + H[31]*P[71];
   const float x145 = x144*x48;
   const float x146 = x144*x55;
   const float x147 = x140*x57;
   const float x148 = x142*x59;
   const float x149 = x144*x62;
   const float x150 = x142*x64;
   const float x151 = x140*x66;
   const float x152 = H[20]*P[51] + H[21]*P[59] + H[22]*P[66] + H[23]*P[72];
   const float x153 = x152*x42;
   const float x154 = H[24]*P[51] + H[25]*P[59] + H[26]*P[66] + H[27]*P[72];
   const float x155 = x154*x45;
   const float x156 = H[28]*P[51] + H[29]*P[59] + H[30]*P[66] + H[31]*P[72];
   const float x157 = x156*x48;
   const float x158 = x156*x55;
   const float x159 = x152*x57;
   const float x160 = x154*x59;
   const float x161 = x156*x62;
   const float x162 = x154*x64;
   const float x163 = x152*x66;
   const float x164 = H[20]*P[52] + H[21]*P[60] + H[22]*P[67] + H[23]*P[73];
   const float x165 = x164*x42;
   const float x166 = H[24]*P[52] + H[25]*P[60] + H[26]*P[67] + H[27]*P[73];
   const float x167 = x166*x45;
   const float x168 = H[28]*P[52] + H[29]*P[60] + H[30]*P[67] + H[31]*P[73];
   const float x169 = x168*x48;
   const float x170 = x168*x55;
   const float x171 = x164*x57;
   const float x172 = x166*x59;
   const float x173 = x168*x62;
   const float x174 = x166*x64;
   const float x175 = x164*x66;
   const float x176 = H[20]*P[53] + H[21]*P[61] + H[22]*P[68] + H[23]*P[74];
   const float x177 = x176*x42;
   const float x178 = H[24]*P[53] + H[25]*P[61] + H[26]*P[68] + H[27]*P[74];
   const float x179 = x178*x45;
   const float x180 = H[28]*P[53] + H[29]*P[61] + H[30]*P[68] + H[31]*P[74];
   const float x181 = x180*x48;
   const float x182 = x180*x55;
   const float x183 = x176*x57;
   const float x184 = x178*x59;
   const float x185 = x180*x62;
   const float x186 = x178*x64;
   const float x187 = x176*x66;
   const float x188 = H[20]*P[54] + H[21]*P[62] + H[22]*P[69] + H[23]*P[75];
   const float x189 = x188*x42;
   const float x190 = H[24]*P[54] + H[25]*P[62] + H[26]*P[69] + H[27]*P[75];
   const float x191 = x190*x45;
   const float x192 = H[28]*P[54] + H[29]*P[62] + H[30]*P[69] + H[31]*P[75];
   const float x193 = x192*x48;
   const float x194 = x192*x55;
   const float x195 = x188*x57;
   const float x196 = x190*x59;
   const float x197 = x192*x62;
   const float x198 = x190*x64;
   const float x199 = x188*x66;
   const float x200 = x43 + x46 + x49;
   const float x201 = x56 + x58 + x60;
   const float x202 = x63 + x65 + x67;
   const float x203 = -x200*H[28] - x201*H[24] - x202*H[20];
   const float x204 = -x200*H[29] - x201*H[25] - x202*H[21];
   const float x205 = -x200*H[30] - x201*H[26] - x202*H[22];
   const float x206 = -x200*H[31] - x201*H[27] - x202*H[23];
   const float x207 = x69 + x71 + x73;
   const float x208 = x74 + x75 + x76;
   const float x209 = x77 + x78 + x79;
   const float x210 = -x207*H[28] - x208*H[24] - x209*H[20];
   const float x211 = -x207*H[29] - x208*H[25] - x209*H[21];
   const float x212 = -x207*H[30] - x208*H[26] - x209*H[22];
   const float x213 = -x207*H[31] - x208*H[27] - x209*H[23];
   const float x214 = x81 + x83 + x85;
   const float x215 = x86 + x87 + x88;
   const float x216 = x89 + x90 + x91;
   const float x217 = -x214*H[28] - x215*H[24] - x216*H[20];
   const float x218 = -x214*H[29] - x215*H[25] - x216*H[21];
   const float x219 = -x214*H[30] - x215*H[26] - x216*H[22];
   const float x220 = -x214*H[31] - x215*H[27] - x216*H[23];
   const float x221 = x93 + x95 + x97;
   const float x222 = x100 + x98 + x99;
   const float x223 = x101 + x102 + x103;
   const float x224 = -x221*H[28] - x222*H[24] - x223*H[20];
   const float x225 = -x221*H[29] - x222*H[25] - x223*H[21];
   const float x226 = -x221*H[30] - x222*H[26] - x223*H[22];
   const float x227 = -x221*H[31] - x222*H[27] - x223*H[23];
   const float x228 = x104 + x105 + x106;
   const float x229 = x107 + x108 + x109;
   const float x230 = x110 + x111 + x112;
   const float x231 = -x228*H[29] - x229*H[25] - x230*H[21];
   const float x232 = -x228*H[30] - x229*H[26] - x230*H[22];
   const float x233 = -x228*H[31] - x229*H[27] - x230*H[23];
   const float x234 = -x228*H[28] - x229*H[24] - x230*H[20] + 1;
   const float x235 = x113 + x114 + x115;
   const float x236 = x116 + x117 + x118;
   const float x237 = x119 + x120 + x121;
   const float x238 = -x235*H[28] - x236*H[24] - x237*H[20];
   const float x239 = -x235*H[30] - x236*H[26] - x237*H[22];
   const float x240 = -x235*H[31] - x236*H[27] - x237*H[23];
   const float x241 = -x235*H[29] - x236*H[25] - x237*H[21] + 1;
   const float x242 = x122 + x123 + x124;
   const float x243 = x125 + x126 + x127;
   const float x244 = x128 + x129 + x130;
   const float x245 = -x242*H[28] - x243*H[24] - x244*H[20];
   const float x246 = -x242*H[29] - x243*H[25] - x244*H[21];
   const float x247 = -x242*H[31] - x243*H[27] - x244*H[23];
   const float x248 = -x242*H[30] - x243*H[26] - x244*H[22] + 1;
   const float x249 = x131 + x132 + x133;
   const float x250 = x134 + x135 + x136;
   const float x251 = x137 + x138 + x139;
   const float x252 = -x249*H[28] - x250*H[24] - x251*H[20];
   const float x253 = -x249*H[29] - x250*H[25] - x251*H[21];
   const float x254 = -x249*H[30] - x250*H[26] - x251*H[22];
   const float x255 = -x249*H[31] - x250*H[27] - x251*H[23] + 1;
   const float x256 = x141 + x143 + x145;
   const float x257 = x146 + x147 + x148;
   const float x258 = x149 + x150 + x151;
   const float x259 = x153 + x155 + x157;
   const float x260 = x158 + x159 + x160;
   const float x261 = x161 + x162 + x163;
   const float x262 = x165 + x167 + x169;
   const float x263 = x170 + x171 + x172;
   const float x264 = x173 + x174 + x175;
   const float x265 = x177 + x179 + x181;
   const float x266 = x182 + x183 + x184;
   const float x267 = x185 + x186 + x187;
   const float x268 = -x265*H[28] - x266*H[24] - x267*H[20];
   const float x269 = -x265*H[29] - x266*H[25] - x267*H[21];
   const float x270 = -x265*H[30] - x266*H[26] - x267*H[22];
   const float x271 = -x265*H[31] - x266*H[27] - x267*H[23];
   const float x272 = x189 + x191 + x193;
   const float x273 = x194 + x195 + x196;
   const float x274 = x197 + x198 + x199;
   xnew[0] = x12*(-x43 - x46 - x49) + x54*(-x56 - x58 - x60) + x61*(-x63 - x65 - x67) + x[0];
   xnew[1] = x12*(-x69 - x71 - x73) + x54*(-x74 - x75 - x76) + x61*(-x77 - x78 - x79) + x[1];
   xnew[2] = x12*(-x81 - x83 - x85) + x54*(-x86 - x87 - x88) + x61*(-x89 - x90 - x91) + x[2];
   xnew[3] = x12*(-x93 - x95 - x97) + x54*(-x100 - x98 - x99) + x61*(-x101 - x102 - x103) + x[3];
   xnew[4] = x12*(-x104 - x105 - x106) + x54*(-x107 - x108 - x109) + x61*(-x110 - x111 - x112) + x[4];
   xnew[5] = x12*(-x113 - x114 - x115) + x54*(-x116 - x117 - x118) + x61*(-x119 - x120 - x121) + x[5];
   xnew[6] = x12*(-x122 - x123 - x124) + x54*(-x125 - x126 - x127) + x61*(-x128 - x129 - x130) + x[6];
   xnew[7] = x12*(-x131 - x132 - x133) + x54*(-x134 - x135 - x136) + x61*(-x137 - x138 - x139) + x[7];
   xnew[8] = x12*(-x141 - x143 - x145) + x54*(-x146 - x147 - x148) + x61*(-x149 - x150 - x151) + x[8];
   xnew[9] = x12*(-x153 - x155 - x157) + x54*(-x158 - x159 - x160) + x61*(-x161 - x162 - x163) + x[9];
   xnew[10] = x12*(-x165 - x167 - x169) + x54*(-x170 - x171 - x172) + x61*(-x173 - x174 - x175) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x12*(-x177 - x179 - x181) + x54*(-x182 - x183 - x184) + x61*(-x185 - x186 - x187) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = x12*(-x189 - x191 - x193) + x54*(-x194 - x195 - x196) + x61*(-x197 - x198 - x199) + x[18];
   Pnew[0] = x203*P[4] + x204*P[5] + x205*P[6] + x206*P[7] + P[0];
   Pnew[1] = x203*P[16] + x204*P[17] + x205*P[18] + x206*P[19] + P[1];
   Pnew[2] = x203*P[27] + x204*P[28] + x205*P[29] + x206*P[30] + P[2];
   Pnew[3] = x203*P[37] + x204*P[38] + x205*P[39] + x206*P[40] + P[3];
   Pnew[4] = x203*P[46] + x204*P[47] + x205*P[48] + x206*P[49] + P[4];
   Pnew[5] = x203*P[47] + x204*P[55] + x205*P[56] + x206*P[57] + P[5];
   Pnew[6] = x203*P[48] + x204*P[56] + x205*P[63] + x206*P[64] + P[6];
   Pnew[7] = x203*P[49] + x204*P[57] + x205*P[64] + x206*P[70] + P[7];
   Pnew[8] = x203*P[50] + x204*P[58] + x205*P[65] + x206*P[71] + P[8];
   Pnew[9] = x203*P[51] + x204*P[59] + x205*P[66] + x206*P[72] + P[9];
   Pnew[10] = x203*P[52] + x204*P[60] + x205*P[67] + x206*P[73] + P[10];
   Pnew[11] = x203*P[53] + x204*P[61] + x205*P[68] + x206*P[74] + P[11];
   Pnew[12] = x203*P[54] + x204*P[62] + x205*P[69] + x206*P[75] + P[12];
   Pnew[13] = x210*P[16] + x211*P[17] + x212*P[18] + x213*P[19] + P[13];
   Pnew[14] = x210*P[27] + x211*P[28] + x212*P[29] + x213*P[30] + P[14];
   Pnew[15] = x210*P[37] + x211*P[38] + x212*P[39] + x213*P[40] + P[15];
   Pnew[16] = x210*P[46] + x211*P[47] + x212*P[48] + x213*P[49] + P[16];
   Pnew[17] = x210*P[47] + x211*P[55] + x212*P[56] + x213*P[57] + P[17];
   Pnew[18] = x210*P[48] + x211*P[56] + x212*P[63] + x213*P[64] + P[18];
   Pnew[19] = x210*P[49] + x211*P[57] + x212*P[64] + x213*P[70] + P[19];
   Pnew[20] = x210*P[50] + x211*P[58] + x212*P[65] + x213*P[71] + P[20];
   Pnew[21] = x210*P[51] + x211*P[59] + x212*P[66] + x213*P[72] + P[21];
   Pnew[22] = x210*P[52] + x211*P[60] + x212*P[67] + x213*P[73] + P[22];
   Pnew[23] = x210*P[53] + x211*P[61] + x212*P[68] + x213*P[74] + P[23];
   Pnew[24] = x210*P[54] + x211*P[62] + x212*P[69] + x213*P[75] + P[24];
   Pnew[25] = x217*P[27] + x218*P[28] + x219*P[29] + x220*P[30] + P[25];
   Pnew[26] = x217*P[37] + x218*P[38] + x219*P[39] + x220*P[40] + P[26];
   Pnew[27] = x217*P[46] + x218*P[47] + x219*P[48] + x220*P[49] + P[27];
   Pnew[28] = x217*P[47] + x218*P[55] + x219*P[56] + x220*P[57] + P[28];
   Pnew[29] = x217*P[48] + x218*P[56] + x219*P[63] + x220*P[64] + P[29];
   Pnew[30] = x217*P[49] + x218*P[57] + x219*P[64] + x220*P[70] + P[30];
   Pnew[31] = x217*P[50] + x218*P[58] + x219*P[65] + x220*P[71] + P[31];
   Pnew[32] = x217*P[51] + x218*P[59] + x219*P[66] + x220*P[72] + P[32];
   Pnew[33] = x217*P[52] + x218*P[60] + x219*P[67] + x220*P[73] + P[33];
   Pnew[34] = x217*P[53] + x218*P[61] + x219*P[68] + x220*P[74] + P[34];
   Pnew[35] = x217*P[54] + x218*P[62] + x219*P[69] + x220*P[75] + P[35];
   Pnew[36] = x224*P[37] + x225*P[38] + x226*P[39] + x227*P[40] + P[36];
   Pnew[37] = x224*P[46] + x225*P[47] + x226*P[48] + x227*P[49] + P[37];
   Pnew[38] = x224*P[47] + x225*P[55] + x226*P[56] + x227*P[57] + P[38];
   Pnew[39] = x224*P[48] + x225*P[56] + x226*P[63] + x227*P[64] + P[39];
   Pnew[40] = x224*P[49] + x225*P[57] + x226*P[64] + x227*P[70] + P[40];
   Pnew[41] = x224*P[50] + x225*P[58] + x226*P[65] + x227*P[71] + P[41];
   Pnew[42] = x224*P[51] + x225*P[59] + x226*P[66] + x227*P[72] + P[42];
   Pnew[43] = x224*P[52] + x225*P[60] + x226*P[67] + x227*P[73] + P[43];
   Pnew[44] = x224*P[53] + x225*P[61] + x226*P[68] + x227*P[74] + P[44];
   Pnew[45] = x224*P[54] + x225*P[62] + x226*P[69] + x227*P[75] + P[45];
   Pnew[46] = x231*P[47] + x232*P[48] + x233*P[49] + x234*P[46];
   Pnew[47] = x231*P[55] + x232*P[56] + x233*P[57] + x234*P[47];
   Pnew[48] = x231*P[56] + x232*P[63] + x233*P[64] + x234*P[48];
   Pnew[49] = x231*P[57] + x232*P[64] + x233*P[70] + x234*P[49];
   Pnew[50] = x231*P[58] + x232*P[65] + x233*P[71] + x234*P[50];
   Pnew[51] = x231*P[59] + x232*P[66] + x233*P[72] + x234*P[51];
   Pnew[52] = x231*P[60] + x232*P[67] + x233*P[73] + x234*P[52];
   Pnew[53] = x231*P[61] + x232*P[68] + x233*P[74] + x234*P[53];
   Pnew[54] = x231*P[62] + x232*P[69] + x233*P[75] + x234*P[54];
   Pnew[55] = x238*P[47] + x239*P[56] + x240*P[57] + x241*P[55];
   Pnew[56] = x238*P[48] + x239*P[63] + x240*P[64] + x241*P[56];
   Pnew[57] = x238*P[49] + x239*P[64] + x240*P[70] + x241*P[57];
   Pnew[58] = x238*P[50] + x239*P[65] + x240*P[71] + x241*P[58];
   Pnew[59] = x238*P[51] + x239*P[66] + x240*P[72] + x241*P[59];
   Pnew[60] = x238*P[52] + x239*P[67] + x240*P[73] + x241*P[60];
   Pnew[61] = x238*P[53] + x239*P[68] + x240*P[74] + x241*P[61];
   Pnew[62] = x238*P[54] + x239*P[69] + x240*P[75] + x241*P[62];
   Pnew[63] = x245*P[48] + x246*P[56] + x247*P[64] + x248*P[63];
   Pnew[64] = x245*P[49] + x246*P[57] + x247*P[70] + x248*P[64];
   Pnew[65] = x245*P[50] + x246*P[58] + x247*P[71] + x248*P[65];
   Pnew[66] = x245*P[51] + x246*P[59] + x247*P[72] + x248*P[66];
   Pnew[67] = x245*P[52] + x246*P[60] + x247*P[73] + x248*P[67];
   Pnew[68] = x245*P[53] + x246*P[61] + x247*P[74] + x248*P[68];
   Pnew[69] = x245*P[54] + x246*P[62] + x247*P[75] + x248*P[69];
   Pnew[70] = x252*P[49] + x253*P[57] + x254*P[64] + x255*P[70];
   Pnew[71] = x252*P[50] + x253*P[58] + x254*P[65] + x255*P[71];
   Pnew[72] = x252*P[51] + x253*P[59] + x254*P[66] + x255*P[72];
   Pnew[73] = x252*P[52] + x253*P[60] + x254*P[67] + x255*P[73];
   Pnew[74] = x252*P[53] + x253*P[61] + x254*P[68] + x255*P[74];
   Pnew[75] = x252*P[54] + x253*P[62] + x254*P[69] + x255*P[75];
   Pnew[76] = (-x256*H[28] - x257*H[24] - x258*H[20])*P[50] + (-x256*H[29] - x257*H[25] - x258*H[21])*P[58] + (-x256*H[30] - x257*H[26] - x258*H[22])*P[65] + (-x256*H[31] - x257*H[27] - x258*H[23])*P[71] + P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = (-x259*H[28] - x260*H[24] - x261*H[20])*P[51] + (-x259*H[29] - x260*H[25] - x261*H[21])*P[59] + (-x259*H[30] - x260*H[26] - x261*H[22])*P[66] + (-x259*H[31] - x260*H[27] - x261*H[23])*P[72] + P[80];
   Pnew[81] = P[81];
   Pnew[82] = (-x262*H[28] - x263*H[24] - x264*H[20])*P[52] + (-x262*H[29] - x263*H[25] - x264*H[21])*P[60] + (-x262*H[30] - x263*H[26] - x264*H[22])*P[67] + (-x262*H[31] - x263*H[27] - x264*H[23])*P[73] + P[82];
   Pnew[83] = P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
   Pnew[87] = P[87];
   Pnew[88] = P[88];
   Pnew[89] = P[89];
   Pnew[90] = P[90];
   Pnew[91] = x268*P[53] + x269*P[61] + x270*P[68] + x271*P[74] + P[91];
   Pnew[92] = x268*P[54] + x269*P[62] + x270*P[69] + x271*P[75] + P[92];
   Pnew[93] = P[93];
   Pnew[94] = P[94];
   Pnew[95] = P[95];
   Pnew[96] = (-x272*H[28] - x273*H[24] - x274*H[20])*P[54] + (-x272*H[29] - x273*H[25] - x274*H[21])*P[62] + (-x272*H[30] - x273*H[26] - x274*H[22])*P[69] + (-x272*H[31] - x273*H[27] - x274*H[23])*P[75] + P[96];
}

void gyro_correction(const float * restrict x, const float * restrict P, const float * restrict gyro,
	const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = (H[17]*H[17]);
   const float x1 = x0*P[76];
   const float x2 = 1.0F/(x1 + R[4]);
   const float x3 = x2*(-gyro[0] + x[8])*H[17];
   const float x4 = (H[18]*H[18]);
   const float x5 = x4*P[80];
   const float x6 = 1.0F/(x5 + R[5]);
   const float x7 = x6*(-gyro[1] + x[9])*H[18];
   const float x8 = (H[19]*H[19]);
   const float x9 = x8*P[82];
   const float x10 = 1.0F/(x9 + R[6]);
   const float x11 = x10*(-gyro[2] + x[10])*H[19];
   const float x12 = x0*x2;
   const float x13 = x4*x6;
   const float x14 = x10*x8;
   const float x15 = x0*x2*P[8];
   const float x16 = x4*x6*P[9];
   const float x17 = x10*x8*P[10];
   const float x18 = x1*x2;
   const float x19 = x5*x6;
   const float x20 = x10*x9;
   const float x21 = x0*x2*P[20];
   const float x22 = x4*x6*P[21];
   const float x23 = x10*x8*P[22];
   const float x24 = x0*x2*P[31];
   const float x25 = x4*x6*P[32];
   const float x26 = x10*x8*P[33];
   const float x27 = x0*x2*P[41];
   const float x28 = x4*x6*P[42];
   const float x29 = x10*x8*P[43];
   const float x30 = x0*x2*P[50];
   const float x31 = x4*x6*P[51];
   const float x32 = x10*x8*P[52];
   const float x33 = x0*x2*P[58];
   const float x34 = x4*x6*P[59];
   const float x35 = x10*x8*P[60];
   const float x36 = -x18 + 1;
   const float x37 = -x19 + 1;
   const float x38 = -x20 + 1;
   xnew[0] = -x11*P[10] - x3*P[8] - x7*P[9] + x[0];
   xnew[1] = -x11*P[22] - x3*P[20] - x7*P[21] + x[1];
   xnew[2] = -x11*P[33] - x3*P[31] - x7*P[32] + x[2];
   xnew[3] = -x11*P[43] - x3*P[41] - x7*P[42] + x[3];
   xnew[4] = -x11*P[52] - x3*P[50] - x7*P[51] + x[4];
   xnew[5] = -x11*P[60] - x3*P[58] - x7*P[59] + x[5];
   xnew[6] = -x11*P[67] - x3*P[65] - x7*P[66] + x[6];
   xnew[7] = -x11*P[73] - x3*P[71] - x7*P[72] + x[7];
   xnew[8] = -x3*P[76] + x[8];
   xnew[9] = -x7*P[80] + x[9];
   xnew[10] = -x11*P[82] + x[10];
   xnew[11] = -x3*P[77] + x[11];
   xnew[12] = -x3*P[78] + x[12];
   xnew[13] = -x11*P[83] + x[13];
   xnew[14] = x[14];
   xnew[15] = -x3*P[79] + x[15];
   xnew[16] = -x7*P[81] + x[16];
   xnew[17] = -x11*P[84] + x[17];
   xnew[18] = x[18];
   Pnew[0] = -x12*(P[8]*P[8]) - x13*(P[9]*P[9]) - x14*(P[10]*P[10]) + P[0];
   Pnew[1] = -x15*P[20] - x16*P[21] - x17*P[22] + P[1];
   Pnew[2] = -x15*P[31] - x16*P[32] - x17*P[33] + P[2];
   Pnew[3] = -x15*P[41] - x16*P[42] - x17*P[43] + P[3];
   Pnew[4] = -x15*P[50] - x16*P[51] - x17*P[52] + P[4];
   Pnew[5] = -x15*P[58] - x16*P[59] - x17*P[60] + P[5];
   Pnew[6] = -x15*P[65] - x16*P[66] - x17*P[67] + P[6];
   Pnew[7] = -x15*P[71] - x16*P[72] - x17*P[73] + P[7];
   Pnew[8] = -x18*P[8] + P[8];
   Pnew[9] = -x19*P[9] + P[9];
   Pnew[10] = -x20*P[10] + P[10];
   Pnew[11] = P[11];
   Pnew[12] = P[12];
   Pnew[13] = -x12*(P[20]*P[20]) - x13*(P[21]*P[21]) - x14*(P[22]*P[22]) + P[13];
   Pnew[14] = -x21*P[31] - x22*P[32] - x23*P[33] + P[14];
   Pnew[15] = -x21*P[41] - x22*P[42] - x23*P[43] + P[15];
   Pnew[16] = -x21*P[50] - x22*P[51] - x23*P[52] + P[16];
   Pnew[17] = -x21*P[58] - x22*P[59] - x23*P[60] + P[17];
   Pnew[18] = -x21*P[65] - x22*P[66] - x23*P[67] + P[18];
   Pnew[19] = -x21*P[71] - x22*P[72] - x23*P[73] + P[19];
   Pnew[20] = -x18*P[20] + P[20];
   Pnew[21] = -x19*P[21] + P[21];
   Pnew[22] = -x20*P[22] + P[22];
   Pnew[23] = P[23];
   Pnew[24] = P[24];
   Pnew[25] = -x12*(P[31]*P[31]) - x13*(P[32]*P[32]) - x14*(P[33]*P[33]) + P[25];
   Pnew[26] = -x24*P[41] - x25*P[42] - x26*P[43] + P[26];
   Pnew[27] = -x24*P[50] - x25*P[51] - x26*P[52] + P[27];
   Pnew[28] = -x24*P[58] - x25*P[59] - x26*P[60] + P[28];
   Pnew[29] = -x24*P[65] - x25*P[66] - x26*P[67] + P[29];
   Pnew[30] = -x24*P[71] - x25*P[72] - x26*P[73] + P[30];
   Pnew[31] = -x18*P[31] + P[31];
   Pnew[32] = -x19*P[32] + P[32];
   Pnew[33] = -x20*P[33] + P[33];
   Pnew[34] = P[34];
   Pnew[35] = P[35];
   Pnew[36] = -x12*(P[41]*P[41]) - x13*(P[42]*P[42]) - x14*(P[43]*P[43]) + P[36];
   Pnew[37] = -x27*P[50] - x28*P[51] - x29*P[52] + P[37];
   Pnew[38] = -x27*P[58] - x28*P[59] - x29*P[60] + P[38];
   Pnew[39] = -x27*P[65] - x28*P[66] - x29*P[67] + P[39];
   Pnew[40] = -x27*P[71] - x28*P[72] - x29*P[73] + P[40];
   Pnew[41] = -x18*P[41] + P[41];
   Pnew[42] = -x19*P[42] + P[42];
   Pnew[43] = -x20*P[43] + P[43];
   Pnew[44] = P[44];
   Pnew[45] = P[45];
   Pnew[46] = -x12*(P[50]*P[50]) - x13*(P[51]*P[51]) - x14*(P[52]*P[52]) + P[46];
   Pnew[47] = -x30*P[58] - x31*P[59] - x32*P[60] + P[47];
   Pnew[48] = -x30*P[65] - x31*P[66] - x32*P[67] + P[48];
   Pnew[49] = -x30*P[71] - x31*P[72] - x32*P[73] + P[49];
   Pnew[50] = -x18*P[50] + P[50];
   Pnew[51] = -x19*P[51] + P[51];
   Pnew[52] = -x20*P[52] + P[52];
   Pnew[53] = P[53];
   Pnew[54] = P[54];
   Pnew[55] = -x12*(P[58]*P[58]) - x13*(P[59]*P[59]) - x14*(P[60]*P[60]) + P[55];
   Pnew[56] = -x33*P[65] - x34*P[66] - x35*P[67] + P[56];
   Pnew[57] = -x33*P[71] - x34*P[72] - x35*P[73] + P[57];
   Pnew[58] = -x18*P[58] + P[58];
   Pnew[59] = -x19*P[59] + P[59];
   Pnew[60] = -x20*P[60] + P[60];
   Pnew[61] = P[61];
   Pnew[62] = P[62];
   Pnew[63] = -x12*(P[65]*P[65]) - x13*(P[66]*P[66]) - x14*(P[67]*P[67]) + P[63];
   Pnew[64] = -x12*P[65]*P[71] - x13*P[66]*P[72] - x14*P[67]*P[73] + P[64];
   Pnew[65] = -x18*P[65] + P[65];
   Pnew[66] = -x19*P[66] + P[66];
   Pnew[67] = -x20*P[67] + P[67];
   Pnew[68] = P[68];
   Pnew[69] = P[69];
   Pnew[70] = -x12*(P[71]*P[71]) - x13*(P[72]*P[72]) - x14*(P[73]*P[73]) + P[70];
   Pnew[71] = -x18*P[71] + P[71];
   Pnew[72] = -x19*P[72] + P[72];
   Pnew[73] = -x20*P[73] + P[73];
   Pnew[74] = P[74];
   Pnew[75] = P[75];
   Pnew[76] = x36*P[76];
   Pnew[77] = x36*P[77];
   Pnew[78] = x36*P[78];
   Pnew[79] = x36*P[79];
   Pnew[80] = x37*P[80];
   Pnew[81] = x37*P[81];
   Pnew[82] = x38*P[82];
   Pnew[83] = x38*P[83];
   Pnew[84] = x38*P[84];
   Pnew[85] = -x12*(P[77]*P[77]) + P[85];
   Pnew[86] = -x12*P[77]*P[79] + P[86];
   Pnew[87] = -x12*(P[78]*P[78]) + P[87];
   Pnew[88] = P[88];
   Pnew[89] = -x14*(P[83]*P[83]) + P[89];
   Pnew[90] = -x14*P[83]*P[84] + P[90];
   Pnew[91] = P[91];
   Pnew[92] = P[92];
   Pnew[93] = -x12*(P[79]*P[79]) + P[93];
   Pnew[94] = -x13*(P[81]*P[81]) + P[94];
   Pnew[95] = -x14*(P[84]*P[84]) + P[95];
   Pnew[96] = P[96];
}

void accel_correction(const float *restrict x, const float *restrict P, const float *restrict accel,
        const float *restrict H, const float *restrict R, const float *restrict param, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = 2*x[3];
   const float x1 = x[4]*x[7];
   const float x2 = x[5]*x[6];
   const float x3 = (x[4]*x[4]) - (x[7]*x[7]);
   const float x4 = (x[6]*x[6]);
   const float x5 = (x[5]*x[5]);
   const float x6 = -(x0*(x[4]*x[5] + x[6]*x[7]) - 2*(x1 - x2)*x[1] + (x3 + x4 - x5)*x[2])*param[6] - accel[1];
   const float x7 = H[1]*P[1] + H[2]*P[2] + H[3]*P[3] + H[4]*P[4] + H[5]*P[5] + H[6]*P[6] + H[7]*P[7];
   const float x8 = H[1]*P[23] + H[2]*P[34] + H[3]*P[44] + H[4]*P[53] + H[5]*P[61] + H[6]*P[68] + H[7]*P[74];
   const float x9 = H[1]*P[24] + H[2]*P[35] + H[3]*P[45] + H[4]*P[54] + H[5]*P[62] + H[6]*P[69] + H[7]*P[75];
   const float x10 = x8*H[15] + x9*H[16];
   const float x11 = H[15]*P[23];
   const float x12 = H[16]*P[24];
   const float x13 = x11 + x12;
   const float x14 = H[15]*P[34];
   const float x15 = H[16]*P[35];
   const float x16 = x14 + x15;
   const float x17 = H[15]*P[44];
   const float x18 = H[16]*P[45];
   const float x19 = x17 + x18;
   const float x20 = H[15]*P[53];
   const float x21 = H[16]*P[54];
   const float x22 = x20 + x21;
   const float x23 = H[15]*P[61];
   const float x24 = H[16]*P[62];
   const float x25 = x23 + x24;
   const float x26 = H[15]*P[68];
   const float x27 = H[16]*P[69];
   const float x28 = x26 + x27;
   const float x29 = H[15]*P[74];
   const float x30 = H[16]*P[75];
   const float x31 = x29 + x30;
   const float x32 = x13*H[8] + x16*H[9] + x19*H[10] + x22*H[11] + x25*H[12] + x28*H[13] + x31*H[14];
   const float x33 = x10*x32;
   const float x34 = H[15]*P[91];
   const float x35 = H[16]*P[92];
   const float x36 = x34 + x35;
   const float x37 = H[15]*P[92];
   const float x38 = H[16]*P[96];
   const float x39 = x37 + x38;
   const float x40 = x36*H[15] + x39*H[16] + R[3];
   const float x41 = H[1]*P[13] + H[2]*P[14] + H[3]*P[15] + H[4]*P[16] + H[5]*P[17] + H[6]*P[18] + H[7]*P[19];
   const float x42 = H[1]*P[14] + H[2]*P[25] + H[3]*P[26] + H[4]*P[27] + H[5]*P[28] + H[6]*P[29] + H[7]*P[30];
   const float x43 = H[1]*P[15] + H[2]*P[26] + H[3]*P[36] + H[4]*P[37] + H[5]*P[38] + H[6]*P[39] + H[7]*P[40];
   const float x44 = H[1]*P[16] + H[2]*P[27] + H[3]*P[37] + H[4]*P[46] + H[5]*P[47] + H[6]*P[48] + H[7]*P[49];
   const float x45 = H[1]*P[17] + H[2]*P[28] + H[3]*P[38] + H[4]*P[47] + H[5]*P[55] + H[6]*P[56] + H[7]*P[57];
   const float x46 = H[1]*P[18] + H[2]*P[29] + H[3]*P[39] + H[4]*P[48] + H[5]*P[56] + H[6]*P[63] + H[7]*P[64];
   const float x47 = H[1]*P[19] + H[2]*P[30] + H[3]*P[40] + H[4]*P[49] + H[5]*P[57] + H[6]*P[64] + H[7]*P[70];
   const float x48 = x41*H[8] + x42*H[9] + x43*H[10] + x44*H[11] + x45*H[12] + x46*H[13] + x47*H[14];
   const float x49 = x40*x48;
   const float x50 = H[8]*P[13] + H[9]*P[14] + H[10]*P[15] + H[11]*P[16] + H[12]*P[17] + H[13]*P[18] + H[14]*P[19];
   const float x51 = H[8]*P[14] + H[9]*P[25] + H[10]*P[26] + H[11]*P[27] + H[12]*P[28] + H[13]*P[29] + H[14]*P[30];
   const float x52 = H[8]*P[15] + H[9]*P[26] + H[10]*P[36] + H[11]*P[37] + H[12]*P[38] + H[13]*P[39] + H[14]*P[40];
   const float x53 = H[8]*P[16] + H[9]*P[27] + H[10]*P[37] + H[11]*P[46] + H[12]*P[47] + H[13]*P[48] + H[14]*P[49];
   const float x54 = H[8]*P[17] + H[9]*P[28] + H[10]*P[38] + H[11]*P[47] + H[12]*P[55] + H[13]*P[56] + H[14]*P[57];
   const float x55 = H[8]*P[18] + H[9]*P[29] + H[10]*P[39] + H[11]*P[48] + H[12]*P[56] + H[13]*P[63] + H[14]*P[64];
   const float x56 = H[8]*P[19] + H[9]*P[30] + H[10]*P[40] + H[11]*P[49] + H[12]*P[57] + H[13]*P[64] + H[14]*P[70];
   const float x57 = x50*H[1] + x51*H[2] + x52*H[3] + x53*H[4] + x54*H[5] + x55*H[6] + x56*H[7];
   const float x58 = H[8]*P[23] + H[9]*P[34] + H[10]*P[44] + H[11]*P[53] + H[12]*P[61] + H[13]*P[68] + H[14]*P[74];
   const float x59 = H[8]*P[24] + H[9]*P[35] + H[10]*P[45] + H[11]*P[54] + H[12]*P[62] + H[13]*P[69] + H[14]*P[75];
   const float x60 = x58*H[15] + x59*H[16];
   const float x61 = x13*H[1] + x16*H[2] + x19*H[3] + x22*H[4] + x25*H[5] + x28*H[6] + x31*H[7];
   const float x62 = x60*x61;
   const float x63 = x10*x61;
   const float x64 = x50*H[8] + x51*H[9] + x52*H[10] + x53*H[11] + x54*H[12] + x55*H[13] + x56*H[14] + R[2];
   const float x65 = x32*x60;
   const float x66 = x41*H[1] + x42*H[2] + x43*H[3] + x44*H[4] + x45*H[5] + x46*H[6] + x47*H[7] + R[1];
   const float x67 = x40*x66;
   const float x68 = 1.0F/(x33*x57 + x48*x62 - x49*x57 - x63*x64 + x64*x67 - x65*x66);
   const float x69 = x68*(x33 - x49);
   const float x70 = x69*x7;
   const float x71 = H[8]*P[1] + H[9]*P[2] + H[10]*P[3] + H[11]*P[4] + H[12]*P[5] + H[13]*P[6] + H[14]*P[7];
   const float x72 = x68*(-x63 + x67);
   const float x73 = x71*x72;
   const float x74 = H[15]*P[11];
   const float x75 = H[16]*P[12];
   const float x76 = x74 + x75;
   const float x77 = x68*(-x32*x66 + x48*x61);
   const float x78 = x76*x77;
   const float x79 = -(-x0*(x[4]*x[6] - x[5]*x[7]) + 2*(x1 + x2)*x[2] + (x3 - x4 + x5)*x[1])*param[6] - accel[0];
   const float x80 = x68*(-x40*x57 + x62);
   const float x81 = x71*x80;
   const float x82 = x68*(x40*x64 - x65);
   const float x83 = x7*x82;
   const float x84 = x68*(x32*x57 - x61*x64);
   const float x85 = x76*x84;
   const float x86 = -accel[2] - param[0]*x[14]*x[18];
   const float x87 = x68*(x10*x57 - x60*x66);
   const float x88 = x71*x87;
   const float x89 = x68*(-x10*x64 + x48*x60);
   const float x90 = x7*x89;
   const float x91 = x68*(-x48*x57 + x64*x66);
   const float x92 = x76*x91;
   const float x93 = x41*x69;
   const float x94 = x50*x72;
   const float x95 = x13*x77;
   const float x96 = x50*x80;
   const float x97 = x41*x82;
   const float x98 = x13*x84;
   const float x99 = x50*x87;
   const float x100 = x41*x89;
   const float x101 = x13*x91;
   const float x102 = x42*x69;
   const float x103 = x51*x72;
   const float x104 = x16*x77;
   const float x105 = x51*x80;
   const float x106 = x42*x82;
   const float x107 = x16*x84;
   const float x108 = x51*x87;
   const float x109 = x42*x89;
   const float x110 = x16*x91;
   const float x111 = x43*x69;
   const float x112 = x52*x72;
   const float x113 = x19*x77;
   const float x114 = x52*x80;
   const float x115 = x43*x82;
   const float x116 = x19*x84;
   const float x117 = x52*x87;
   const float x118 = x43*x89;
   const float x119 = x19*x91;
   const float x120 = x44*x69;
   const float x121 = x53*x72;
   const float x122 = x22*x77;
   const float x123 = x53*x80;
   const float x124 = x44*x82;
   const float x125 = x22*x84;
   const float x126 = x53*x87;
   const float x127 = x44*x89;
   const float x128 = x22*x91;
   const float x129 = x45*x69;
   const float x130 = x54*x72;
   const float x131 = x25*x77;
   const float x132 = x54*x80;
   const float x133 = x45*x82;
   const float x134 = x25*x84;
   const float x135 = x54*x87;
   const float x136 = x45*x89;
   const float x137 = x25*x91;
   const float x138 = x46*x69;
   const float x139 = x55*x72;
   const float x140 = x28*x77;
   const float x141 = x55*x80;
   const float x142 = x46*x82;
   const float x143 = x28*x84;
   const float x144 = x55*x87;
   const float x145 = x46*x89;
   const float x146 = x28*x91;
   const float x147 = x47*x69;
   const float x148 = x56*x72;
   const float x149 = x31*x77;
   const float x150 = x56*x80;
   const float x151 = x47*x82;
   const float x152 = x31*x84;
   const float x153 = x56*x87;
   const float x154 = x47*x89;
   const float x155 = x31*x91;
   const float x156 = H[1]*P[20] + H[2]*P[31] + H[3]*P[41] + H[4]*P[50] + H[5]*P[58] + H[6]*P[65] + H[7]*P[71];
   const float x157 = x156*x69;
   const float x158 = H[8]*P[20] + H[9]*P[31] + H[10]*P[41] + H[11]*P[50] + H[12]*P[58] + H[13]*P[65] + H[14]*P[71];
   const float x159 = x158*x72;
   const float x160 = x158*x80;
   const float x161 = x156*x82;
   const float x162 = H[1]*P[21] + H[2]*P[32] + H[3]*P[42] + H[4]*P[51] + H[5]*P[59] + H[6]*P[66] + H[7]*P[72];
   const float x163 = x162*x69;
   const float x164 = H[8]*P[21] + H[9]*P[32] + H[10]*P[42] + H[11]*P[51] + H[12]*P[59] + H[13]*P[66] + H[14]*P[72];
   const float x165 = x164*x72;
   const float x166 = x164*x80;
   const float x167 = x162*x82;
   const float x168 = H[1]*P[22] + H[2]*P[33] + H[3]*P[43] + H[4]*P[52] + H[5]*P[60] + H[6]*P[67] + H[7]*P[73];
   const float x169 = x168*x69;
   const float x170 = H[8]*P[22] + H[9]*P[33] + H[10]*P[43] + H[11]*P[52] + H[12]*P[60] + H[13]*P[67] + H[14]*P[73];
   const float x171 = x170*x72;
   const float x172 = x170*x80;
   const float x173 = x168*x82;
   const float x174 = x69*x8;
   const float x175 = x58*x72;
   const float x176 = x36*x77;
   const float x177 = x58*x80;
   const float x178 = x8*x82;
   const float x179 = x36*x84;
   const float x180 = x58*x87;
   const float x181 = x8*x89;
   const float x182 = x36*x91;
   const float x183 = x69*x9;
   const float x184 = x59*x72;
   const float x185 = x39*x77;
   const float x186 = x59*x80;
   const float x187 = x82*x9;
   const float x188 = x39*x84;
   const float x189 = x59*x87;
   const float x190 = x89*x9;
   const float x191 = x39*x91;
   const float x192 = x88 + x90 + x92;
   const float x193 = x70 + x73 + x78;
   const float x194 = x81 + x83 + x85;
   const float x195 = -x193*H[8] - x194*H[1];
   const float x196 = -x193*H[9] - x194*H[2];
   const float x197 = -x193*H[10] - x194*H[3];
   const float x198 = -x193*H[11] - x194*H[4];
   const float x199 = -x193*H[12] - x194*H[5];
   const float x200 = -x193*H[13] - x194*H[6];
   const float x201 = -x193*H[14] - x194*H[7];
   const float x202 = x100 + x101 + x99;
   const float x203 = x93 + x94 + x95;
   const float x204 = x96 + x97 + x98;
   const float x205 = -x203*H[9] - x204*H[2];
   const float x206 = -x203*H[10] - x204*H[3];
   const float x207 = -x203*H[11] - x204*H[4];
   const float x208 = -x203*H[12] - x204*H[5];
   const float x209 = -x203*H[13] - x204*H[6];
   const float x210 = -x203*H[14] - x204*H[7];
   const float x211 = -x203*H[8] - x204*H[1] + 1;
   const float x212 = x108 + x109 + x110;
   const float x213 = x102 + x103 + x104;
   const float x214 = x105 + x106 + x107;
   const float x215 = -x213*H[8] - x214*H[1];
   const float x216 = -x213*H[10] - x214*H[3];
   const float x217 = -x213*H[11] - x214*H[4];
   const float x218 = -x213*H[12] - x214*H[5];
   const float x219 = -x213*H[13] - x214*H[6];
   const float x220 = -x213*H[14] - x214*H[7];
   const float x221 = -x213*H[9] - x214*H[2] + 1;
   const float x222 = x117 + x118 + x119;
   const float x223 = x111 + x112 + x113;
   const float x224 = x114 + x115 + x116;
   const float x225 = -x223*H[8] - x224*H[1];
   const float x226 = -x223*H[9] - x224*H[2];
   const float x227 = -x223*H[11] - x224*H[4];
   const float x228 = -x223*H[12] - x224*H[5];
   const float x229 = -x223*H[13] - x224*H[6];
   const float x230 = -x223*H[14] - x224*H[7];
   const float x231 = -x223*H[10] - x224*H[3] + 1;
   const float x232 = x126 + x127 + x128;
   const float x233 = x120 + x121 + x122;
   const float x234 = x123 + x124 + x125;
   const float x235 = -x233*H[8] - x234*H[1];
   const float x236 = -x233*H[9] - x234*H[2];
   const float x237 = -x233*H[10] - x234*H[3];
   const float x238 = -x233*H[12] - x234*H[5];
   const float x239 = -x233*H[13] - x234*H[6];
   const float x240 = -x233*H[14] - x234*H[7];
   const float x241 = -x233*H[11] - x234*H[4] + 1;
   const float x242 = x135 + x136 + x137;
   const float x243 = x129 + x130 + x131;
   const float x244 = x132 + x133 + x134;
   const float x245 = -x243*H[8] - x244*H[1];
   const float x246 = -x243*H[9] - x244*H[2];
   const float x247 = -x243*H[10] - x244*H[3];
   const float x248 = -x243*H[11] - x244*H[4];
   const float x249 = -x243*H[13] - x244*H[6];
   const float x250 = -x243*H[14] - x244*H[7];
   const float x251 = -x243*H[12] - x244*H[5] + 1;
   const float x252 = x144 + x145 + x146;
   const float x253 = x138 + x139 + x140;
   const float x254 = x141 + x142 + x143;
   const float x255 = -x253*H[8] - x254*H[1];
   const float x256 = -x253*H[9] - x254*H[2];
   const float x257 = -x253*H[10] - x254*H[3];
   const float x258 = -x253*H[11] - x254*H[4];
   const float x259 = -x253*H[12] - x254*H[5];
   const float x260 = -x253*H[14] - x254*H[7];
   const float x261 = -x253*H[13] - x254*H[6] + 1;
   const float x262 = x153 + x154 + x155;
   const float x263 = x147 + x148 + x149;
   const float x264 = x150 + x151 + x152;
   const float x265 = -x263*H[8] - x264*H[1];
   const float x266 = -x263*H[9] - x264*H[2];
   const float x267 = -x263*H[10] - x264*H[3];
   const float x268 = -x263*H[11] - x264*H[4];
   const float x269 = -x263*H[12] - x264*H[5];
   const float x270 = -x263*H[13] - x264*H[6];
   const float x271 = -x263*H[14] - x264*H[7] + 1;
   const float x272 = x157 + x159;
   const float x273 = x160 + x161;
   const float x274 = x163 + x165;
   const float x275 = x166 + x167;
   const float x276 = x169 + x171;
   const float x277 = x172 + x173;
   const float x278 = x180 + x181 + x182;
   const float x279 = -x278*H[15] + 1;
   const float x280 = x174 + x175 + x176;
   const float x281 = x177 + x178 + x179;
   const float x282 = -x280*H[8] - x281*H[1];
   const float x283 = -x280*H[9] - x281*H[2];
   const float x284 = -x280*H[10] - x281*H[3];
   const float x285 = -x280*H[11] - x281*H[4];
   const float x286 = -x280*H[12] - x281*H[5];
   const float x287 = -x280*H[13] - x281*H[6];
   const float x288 = -x280*H[14] - x281*H[7];
   const float x289 = x189 + x190 + x191;
   const float x290 = x183 + x184 + x185;
   const float x291 = x186 + x187 + x188;
   xnew[0] = x6*(-x70 - x73 - x78) + x79*(-x81 - x83 - x85) + x86*(-x88 - x90 - x92) + x[0];
   xnew[1] = x6*(-x93 - x94 - x95) + x79*(-x96 - x97 - x98) + x86*(-x100 - x101 - x99) + x[1];
   xnew[2] = x6*(-x102 - x103 - x104) + x79*(-x105 - x106 - x107) + x86*(-x108 - x109 - x110) + x[2];
   xnew[3] = x6*(-x111 - x112 - x113) + x79*(-x114 - x115 - x116) + x86*(-x117 - x118 - x119) + x[3];
   xnew[4] = x6*(-x120 - x121 - x122) + x79*(-x123 - x124 - x125) + x86*(-x126 - x127 - x128) + x[4];
   xnew[5] = x6*(-x129 - x130 - x131) + x79*(-x132 - x133 - x134) + x86*(-x135 - x136 - x137) + x[5];
   xnew[6] = x6*(-x138 - x139 - x140) + x79*(-x141 - x142 - x143) + x86*(-x144 - x145 - x146) + x[6];
   xnew[7] = x6*(-x147 - x148 - x149) + x79*(-x150 - x151 - x152) + x86*(-x153 - x154 - x155) + x[7];
   xnew[8] = x6*(-x157 - x159) + x79*(-x160 - x161) + x86*(-x156*x89 - x158*x87) + x[8];
   xnew[9] = x6*(-x163 - x165) + x79*(-x166 - x167) + x86*(-x162*x89 - x164*x87) + x[9];
   xnew[10] = x6*(-x169 - x171) + x79*(-x172 - x173) + x86*(-x168*x89 - x170*x87) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x6*(-x174 - x175 - x176) + x79*(-x177 - x178 - x179) + x86*(-x180 - x181 - x182) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = x6*(-x183 - x184 - x185) + x79*(-x186 - x187 - x188) + x86*(-x189 - x190 - x191) + x[18];
   Pnew[0] = -x192*x74 - x192*x75 + x195*P[1] + x196*P[2] + x197*P[3] + x198*P[4] + x199*P[5] + x200*P[6] + x201*P[7] + P[0];
   Pnew[1] = -x11*x192 - x12*x192 + x195*P[13] + x196*P[14] + x197*P[15] + x198*P[16] + x199*P[17] + x200*P[18] + x201*P[19] + P[1];
   Pnew[2] = -x14*x192 - x15*x192 + x195*P[14] + x196*P[25] + x197*P[26] + x198*P[27] + x199*P[28] + x200*P[29] + x201*P[30] + P[2];
   Pnew[3] = -x17*x192 - x18*x192 + x195*P[15] + x196*P[26] + x197*P[36] + x198*P[37] + x199*P[38] + x200*P[39] + x201*P[40] + P[3];
   Pnew[4] = -x192*x20 - x192*x21 + x195*P[16] + x196*P[27] + x197*P[37] + x198*P[46] + x199*P[47] + x200*P[48] + x201*P[49] + P[4];
   Pnew[5] = -x192*x23 - x192*x24 + x195*P[17] + x196*P[28] + x197*P[38] + x198*P[47] + x199*P[55] + x200*P[56] + x201*P[57] + P[5];
   Pnew[6] = -x192*x26 - x192*x27 + x195*P[18] + x196*P[29] + x197*P[39] + x198*P[48] + x199*P[56] + x200*P[63] + x201*P[64] + P[6];
   Pnew[7] = -x192*x29 - x192*x30 + x195*P[19] + x196*P[30] + x197*P[40] + x198*P[49] + x199*P[57] + x200*P[64] + x201*P[70] + P[7];
   Pnew[8] = x195*P[20] + x196*P[31] + x197*P[41] + x198*P[50] + x199*P[58] + x200*P[65] + x201*P[71] + P[8];
   Pnew[9] = x195*P[21] + x196*P[32] + x197*P[42] + x198*P[51] + x199*P[59] + x200*P[66] + x201*P[72] + P[9];
   Pnew[10] = x195*P[22] + x196*P[33] + x197*P[43] + x198*P[52] + x199*P[60] + x200*P[67] + x201*P[73] + P[10];
   Pnew[11] = -x192*x34 - x192*x35 + x195*P[23] + x196*P[34] + x197*P[44] + x198*P[53] + x199*P[61] + x200*P[68] + x201*P[74] + P[11];
   Pnew[12] = -x192*x37 - x192*x38 + x195*P[24] + x196*P[35] + x197*P[45] + x198*P[54] + x199*P[62] + x200*P[69] + x201*P[75] + P[12];
   Pnew[13] = -x11*x202 - x12*x202 + x205*P[14] + x206*P[15] + x207*P[16] + x208*P[17] + x209*P[18] + x210*P[19] + x211*P[13];
   Pnew[14] = -x14*x202 - x15*x202 + x205*P[25] + x206*P[26] + x207*P[27] + x208*P[28] + x209*P[29] + x210*P[30] + x211*P[14];
   Pnew[15] = -x17*x202 - x18*x202 + x205*P[26] + x206*P[36] + x207*P[37] + x208*P[38] + x209*P[39] + x210*P[40] + x211*P[15];
   Pnew[16] = -x20*x202 - x202*x21 + x205*P[27] + x206*P[37] + x207*P[46] + x208*P[47] + x209*P[48] + x210*P[49] + x211*P[16];
   Pnew[17] = -x202*x23 - x202*x24 + x205*P[28] + x206*P[38] + x207*P[47] + x208*P[55] + x209*P[56] + x210*P[57] + x211*P[17];
   Pnew[18] = -x202*x26 - x202*x27 + x205*P[29] + x206*P[39] + x207*P[48] + x208*P[56] + x209*P[63] + x210*P[64] + x211*P[18];
   Pnew[19] = -x202*x29 - x202*x30 + x205*P[30] + x206*P[40] + x207*P[49] + x208*P[57] + x209*P[64] + x210*P[70] + x211*P[19];
   Pnew[20] = x205*P[31] + x206*P[41] + x207*P[50] + x208*P[58] + x209*P[65] + x210*P[71] + x211*P[20];
   Pnew[21] = x205*P[32] + x206*P[42] + x207*P[51] + x208*P[59] + x209*P[66] + x210*P[72] + x211*P[21];
   Pnew[22] = x205*P[33] + x206*P[43] + x207*P[52] + x208*P[60] + x209*P[67] + x210*P[73] + x211*P[22];
   Pnew[23] = -x202*x34 - x202*x35 + x205*P[34] + x206*P[44] + x207*P[53] + x208*P[61] + x209*P[68] + x210*P[74] + x211*P[23];
   Pnew[24] = -x202*x37 - x202*x38 + x205*P[35] + x206*P[45] + x207*P[54] + x208*P[62] + x209*P[69] + x210*P[75] + x211*P[24];
   Pnew[25] = -x14*x212 - x15*x212 + x215*P[14] + x216*P[26] + x217*P[27] + x218*P[28] + x219*P[29] + x220*P[30] + x221*P[25];
   Pnew[26] = -x17*x212 - x18*x212 + x215*P[15] + x216*P[36] + x217*P[37] + x218*P[38] + x219*P[39] + x220*P[40] + x221*P[26];
   Pnew[27] = -x20*x212 - x21*x212 + x215*P[16] + x216*P[37] + x217*P[46] + x218*P[47] + x219*P[48] + x220*P[49] + x221*P[27];
   Pnew[28] = -x212*x23 - x212*x24 + x215*P[17] + x216*P[38] + x217*P[47] + x218*P[55] + x219*P[56] + x220*P[57] + x221*P[28];
   Pnew[29] = -x212*x26 - x212*x27 + x215*P[18] + x216*P[39] + x217*P[48] + x218*P[56] + x219*P[63] + x220*P[64] + x221*P[29];
   Pnew[30] = -x212*x29 - x212*x30 + x215*P[19] + x216*P[40] + x217*P[49] + x218*P[57] + x219*P[64] + x220*P[70] + x221*P[30];
   Pnew[31] = x215*P[20] + x216*P[41] + x217*P[50] + x218*P[58] + x219*P[65] + x220*P[71] + x221*P[31];
   Pnew[32] = x215*P[21] + x216*P[42] + x217*P[51] + x218*P[59] + x219*P[66] + x220*P[72] + x221*P[32];
   Pnew[33] = x215*P[22] + x216*P[43] + x217*P[52] + x218*P[60] + x219*P[67] + x220*P[73] + x221*P[33];
   Pnew[34] = -x212*x34 - x212*x35 + x215*P[23] + x216*P[44] + x217*P[53] + x218*P[61] + x219*P[68] + x220*P[74] + x221*P[34];
   Pnew[35] = -x212*x37 - x212*x38 + x215*P[24] + x216*P[45] + x217*P[54] + x218*P[62] + x219*P[69] + x220*P[75] + x221*P[35];
   Pnew[36] = -x17*x222 - x18*x222 + x225*P[15] + x226*P[26] + x227*P[37] + x228*P[38] + x229*P[39] + x230*P[40] + x231*P[36];
   Pnew[37] = -x20*x222 - x21*x222 + x225*P[16] + x226*P[27] + x227*P[46] + x228*P[47] + x229*P[48] + x230*P[49] + x231*P[37];
   Pnew[38] = -x222*x23 - x222*x24 + x225*P[17] + x226*P[28] + x227*P[47] + x228*P[55] + x229*P[56] + x230*P[57] + x231*P[38];
   Pnew[39] = -x222*x26 - x222*x27 + x225*P[18] + x226*P[29] + x227*P[48] + x228*P[56] + x229*P[63] + x230*P[64] + x231*P[39];
   Pnew[40] = -x222*x29 - x222*x30 + x225*P[19] + x226*P[30] + x227*P[49] + x228*P[57] + x229*P[64] + x230*P[70] + x231*P[40];
   Pnew[41] = x225*P[20] + x226*P[31] + x227*P[50] + x228*P[58] + x229*P[65] + x230*P[71] + x231*P[41];
   Pnew[42] = x225*P[21] + x226*P[32] + x227*P[51] + x228*P[59] + x229*P[66] + x230*P[72] + x231*P[42];
   Pnew[43] = x225*P[22] + x226*P[33] + x227*P[52] + x228*P[60] + x229*P[67] + x230*P[73] + x231*P[43];
   Pnew[44] = -x222*x34 - x222*x35 + x225*P[23] + x226*P[34] + x227*P[53] + x228*P[61] + x229*P[68] + x230*P[74] + x231*P[44];
   Pnew[45] = -x222*x37 - x222*x38 + x225*P[24] + x226*P[35] + x227*P[54] + x228*P[62] + x229*P[69] + x230*P[75] + x231*P[45];
   Pnew[46] = -x20*x232 - x21*x232 + x235*P[16] + x236*P[27] + x237*P[37] + x238*P[47] + x239*P[48] + x240*P[49] + x241*P[46];
   Pnew[47] = -x23*x232 - x232*x24 + x235*P[17] + x236*P[28] + x237*P[38] + x238*P[55] + x239*P[56] + x240*P[57] + x241*P[47];
   Pnew[48] = -x232*x26 - x232*x27 + x235*P[18] + x236*P[29] + x237*P[39] + x238*P[56] + x239*P[63] + x240*P[64] + x241*P[48];
   Pnew[49] = -x232*x29 - x232*x30 + x235*P[19] + x236*P[30] + x237*P[40] + x238*P[57] + x239*P[64] + x240*P[70] + x241*P[49];
   Pnew[50] = x235*P[20] + x236*P[31] + x237*P[41] + x238*P[58] + x239*P[65] + x240*P[71] + x241*P[50];
   Pnew[51] = x235*P[21] + x236*P[32] + x237*P[42] + x238*P[59] + x239*P[66] + x240*P[72] + x241*P[51];
   Pnew[52] = x235*P[22] + x236*P[33] + x237*P[43] + x238*P[60] + x239*P[67] + x240*P[73] + x241*P[52];
   Pnew[53] = -x232*x34 - x232*x35 + x235*P[23] + x236*P[34] + x237*P[44] + x238*P[61] + x239*P[68] + x240*P[74] + x241*P[53];
   Pnew[54] = -x232*x37 - x232*x38 + x235*P[24] + x236*P[35] + x237*P[45] + x238*P[62] + x239*P[69] + x240*P[75] + x241*P[54];
   Pnew[55] = -x23*x242 - x24*x242 + x245*P[17] + x246*P[28] + x247*P[38] + x248*P[47] + x249*P[56] + x250*P[57] + x251*P[55];
   Pnew[56] = -x242*x26 - x242*x27 + x245*P[18] + x246*P[29] + x247*P[39] + x248*P[48] + x249*P[63] + x250*P[64] + x251*P[56];
   Pnew[57] = -x242*x29 - x242*x30 + x245*P[19] + x246*P[30] + x247*P[40] + x248*P[49] + x249*P[64] + x250*P[70] + x251*P[57];
   Pnew[58] = x245*P[20] + x246*P[31] + x247*P[41] + x248*P[50] + x249*P[65] + x250*P[71] + x251*P[58];
   Pnew[59] = x245*P[21] + x246*P[32] + x247*P[42] + x248*P[51] + x249*P[66] + x250*P[72] + x251*P[59];
   Pnew[60] = x245*P[22] + x246*P[33] + x247*P[43] + x248*P[52] + x249*P[67] + x250*P[73] + x251*P[60];
   Pnew[61] = -x242*x34 - x242*x35 + x245*P[23] + x246*P[34] + x247*P[44] + x248*P[53] + x249*P[68] + x250*P[74] + x251*P[61];
   Pnew[62] = -x242*x37 - x242*x38 + x245*P[24] + x246*P[35] + x247*P[45] + x248*P[54] + x249*P[69] + x250*P[75] + x251*P[62];
   Pnew[63] = -x252*x26 - x252*x27 + x255*P[18] + x256*P[29] + x257*P[39] + x258*P[48] + x259*P[56] + x260*P[64] + x261*P[63];
   Pnew[64] = -x252*x29 - x252*x30 + x255*P[19] + x256*P[30] + x257*P[40] + x258*P[49] + x259*P[57] + x260*P[70] + x261*P[64];
   Pnew[65] = x255*P[20] + x256*P[31] + x257*P[41] + x258*P[50] + x259*P[58] + x260*P[71] + x261*P[65];
   Pnew[66] = x255*P[21] + x256*P[32] + x257*P[42] + x258*P[51] + x259*P[59] + x260*P[72] + x261*P[66];
   Pnew[67] = x255*P[22] + x256*P[33] + x257*P[43] + x258*P[52] + x259*P[60] + x260*P[73] + x261*P[67];
   Pnew[68] = -x252*x34 - x252*x35 + x255*P[23] + x256*P[34] + x257*P[44] + x258*P[53] + x259*P[61] + x260*P[74] + x261*P[68];
   Pnew[69] = -x252*x37 - x252*x38 + x255*P[24] + x256*P[35] + x257*P[45] + x258*P[54] + x259*P[62] + x260*P[75] + x261*P[69];
   Pnew[70] = -x262*x29 - x262*x30 + x265*P[19] + x266*P[30] + x267*P[40] + x268*P[49] + x269*P[57] + x270*P[64] + x271*P[70];
   Pnew[71] = x265*P[20] + x266*P[31] + x267*P[41] + x268*P[50] + x269*P[58] + x270*P[65] + x271*P[71];
   Pnew[72] = x265*P[21] + x266*P[32] + x267*P[42] + x268*P[51] + x269*P[59] + x270*P[66] + x271*P[72];
   Pnew[73] = x265*P[22] + x266*P[33] + x267*P[43] + x268*P[52] + x269*P[60] + x270*P[67] + x271*P[73];
   Pnew[74] = -x262*x34 - x262*x35 + x265*P[23] + x266*P[34] + x267*P[44] + x268*P[53] + x269*P[61] + x270*P[68] + x271*P[74];
   Pnew[75] = -x262*x37 - x262*x38 + x265*P[24] + x266*P[35] + x267*P[45] + x268*P[54] + x269*P[62] + x270*P[69] + x271*P[75];
   Pnew[76] = (-x272*H[8] - x273*H[1])*P[20] + (-x272*H[9] - x273*H[2])*P[31] + (-x272*H[10] - x273*H[3])*P[41] + (-x272*H[11] - x273*H[4])*P[50] + (-x272*H[12] - x273*H[5])*P[58] + (-x272*H[13] - x273*H[6])*P[65] + (-x272*H[14] - x273*H[7])*P[71] + P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = (-x274*H[8] - x275*H[1])*P[21] + (-x274*H[9] - x275*H[2])*P[32] + (-x274*H[10] - x275*H[3])*P[42] + (-x274*H[11] - x275*H[4])*P[51] + (-x274*H[12] - x275*H[5])*P[59] + (-x274*H[13] - x275*H[6])*P[66] + (-x274*H[14] - x275*H[7])*P[72] + P[80];
   Pnew[81] = P[81];
   Pnew[82] = (-x276*H[8] - x277*H[1])*P[22] + (-x276*H[9] - x277*H[2])*P[33] + (-x276*H[10] - x277*H[3])*P[43] + (-x276*H[11] - x277*H[4])*P[52] + (-x276*H[12] - x277*H[5])*P[60] + (-x276*H[13] - x277*H[6])*P[67] + (-x276*H[14] - x277*H[7])*P[73] + P[82];
   Pnew[83] = P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
   Pnew[87] = P[87];
   Pnew[88] = P[88];
   Pnew[89] = P[89];
   Pnew[90] = P[90];
   Pnew[91] = -x278*x35 + x279*P[91] + x282*P[23] + x283*P[34] + x284*P[44] + x285*P[53] + x286*P[61] + x287*P[68] + x288*P[74];
   Pnew[92] = -x278*x38 + x279*P[92] + x282*P[24] + x283*P[35] + x284*P[45] + x285*P[54] + x286*P[62] + x287*P[69] + x288*P[75];
   Pnew[93] = P[93];
   Pnew[94] = P[94];
   Pnew[95] = P[95];
   Pnew[96] = -x289*x37 + (-x289*H[16] + 1)*P[96] + (-x290*H[8] - x291*H[1])*P[24] + (-x290*H[9] - x291*H[2])*P[35] + (-x290*H[10] - x291*H[3])*P[45] + (-x290*H[11] - x291*H[4])*P[54] + (-x290*H[12] - x291*H[5])*P[62] + (-x290*H[13] - x291*H[6])*P[69] + (-x290*H[14] - x291*H[7])*P[75];
}

void covariance_init(float *Pnew) {
	float P0[NUMX] = {1.0e7f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		              1.0e-9f, 1.0e-9f, 1.0e-9f, 1.0e-9f};
     Pnew[0] = P0[0];
   Pnew[1] = 0;
   Pnew[2] = 0;
   Pnew[3] = 0;
   Pnew[4] = 0;
   Pnew[5] = 0;
   Pnew[6] = 0;
   Pnew[7] = 0;
   Pnew[8] = 0;
   Pnew[9] = 0;
   Pnew[10] = 0;
   Pnew[11] = 0;
   Pnew[12] = 0;
   Pnew[13] = P0[1];
   Pnew[14] = 0;
   Pnew[15] = 0;
   Pnew[16] = 0;
   Pnew[17] = 0;
   Pnew[18] = 0;
   Pnew[19] = 0;
   Pnew[20] = 0;
   Pnew[21] = 0;
   Pnew[22] = 0;
   Pnew[23] = 0;
   Pnew[24] = 0;
   Pnew[25] = P0[2];
   Pnew[26] = 0;
   Pnew[27] = 0;
   Pnew[28] = 0;
   Pnew[29] = 0;
   Pnew[30] = 0;
   Pnew[31] = 0;
   Pnew[32] = 0;
   Pnew[33] = 0;
   Pnew[34] = 0;
   Pnew[35] = 0;
   Pnew[36] = P0[3];
   Pnew[37] = 0;
   Pnew[38] = 0;
   Pnew[39] = 0;
   Pnew[40] = 0;
   Pnew[41] = 0;
   Pnew[42] = 0;
   Pnew[43] = 0;
   Pnew[44] = 0;
   Pnew[45] = 0;
   Pnew[46] = P0[4];
   Pnew[47] = 0;
   Pnew[48] = 0;
   Pnew[49] = 0;
   Pnew[50] = 0;
   Pnew[51] = 0;
   Pnew[52] = 0;
   Pnew[53] = 0;
   Pnew[54] = 0;
   Pnew[55] = P0[5];
   Pnew[56] = 0;
   Pnew[57] = 0;
   Pnew[58] = 0;
   Pnew[59] = 0;
   Pnew[60] = 0;
   Pnew[61] = 0;
   Pnew[62] = 0;
   Pnew[63] = P0[6];
   Pnew[64] = 0;
   Pnew[65] = 0;
   Pnew[66] = 0;
   Pnew[67] = 0;
   Pnew[68] = 0;
   Pnew[69] = 0;
   Pnew[70] = P0[7];
   Pnew[71] = 0;
   Pnew[72] = 0;
   Pnew[73] = 0;
   Pnew[74] = 0;
   Pnew[75] = 0;
   Pnew[76] = P0[8];
   Pnew[77] = 0;
   Pnew[78] = 0;
   Pnew[79] = 0;
   Pnew[80] = P0[9];
   Pnew[81] = 0;
   Pnew[82] = P0[10];
   Pnew[83] = 0;
   Pnew[84] = 0;
   Pnew[85] = P0[11];
   Pnew[86] = 0;
   Pnew[87] = P0[12];
   Pnew[88] = 0;
   Pnew[89] = P0[13];
   Pnew[90] = 0;
   Pnew[91] = P0[14];
   Pnew[92] = 0;
   Pnew[93] = P0[15];
   Pnew[94] = P0[16];
   Pnew[95] = P0[17];
   Pnew[96] = P0[18];
}

/**
 * @}
 * @}
 */
