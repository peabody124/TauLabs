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
#define NUMX 23  // number of state variables
#define NUMF 66  // number of non-zero variables stored in the F matrix
#define NUMH 37  // number of non-zero variables stored in the H matrix
#define NUMP 144 // number of non-zero variables in the upper triangular of covariance
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
		float Be[3];
	} __attribute__((__packed__)) params;

	struct {
		float mu;
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

void qcins_unarmed_state(struct qcins_state *qcins_state);

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
	Q[15] = Q[16] = Q[17] = 1e-5f;      // Input bias states
	Q[18] = Q[19] = Q[20] = 1e-5f;      // Output bias states
	Q[21] = 1e-5f;                      // Thrust gain
	Q[22] = 1e-5f;                      // Mu

	qcins_state->armed = false;

	// Store the observation noises
	qcins_state->R[0] = 100;
	qcins_state->R[1] = qcins_state->R[2] = qcins_state->R[3] = 1e6;
	qcins_state->R[4] = qcins_state->R[5] = qcins_state->R[6] = 1e3;
	qcins_state->R[7] = qcins_state->R[8] = qcins_state->R[9] = 1e3;

	qcins_state->init.bias[0] = 0.0f;
	qcins_state->init.bias[1] = 0.0f;
	qcins_state->init.bias[2] = 0.0f;
	qcins_state->init.beta_t = logf(2.0f); // thrust to weight ratio 2:1
	qcins_state->init.mu = -2.0f;

	// Defaults for the parameters
	qcins_state->params.g = 9.81f;
	qcins_state->params.beta_r = 10000.0f * DEG2RAD;
	qcins_state->params.beta_p = 10000.0f * DEG2RAD;
	qcins_state->params.beta_y1 = 1000.0f * DEG2RAD;
	qcins_state->params.beta_y2 = 1000.0f * DEG2RAD;
	qcins_state->params.tau = 0.050f;
	qcins_state->params.Be[0] = 0.36;
	qcins_state->params.Be[1] = -0.03;
	qcins_state->params.Be[2] = 0.93;

	// Initialize some of the states
	qcins_unarmed_state(qcins_state);
	
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

bool qcins_set_init_mu(uintptr_t qcins_handle, const float mu_new)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->init.mu = mu_new;

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

void qcins_unarmed_state(struct qcins_state *qcins_state)
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
	qcins_state->x[21] = qcins_state->init.beta_t; // thrust to weight ratio
	qcins_state->x[22] = qcins_state->init.mu;
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
	const float mag_len = sqrtf(mag[0]*mag[0] + mag[1]*mag[1] + mag[2]*mag[2]) + 1e-10f;
	const float mag_norm[3] = {mag[0] / mag_len, mag[1] / mag_len, mag[2] / mag_len};

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

bool qcins_get_output_bias(uintptr_t qcins_handle, float out_bias[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	out_bias[0] = qcins_state->x[18];
	out_bias[1] = qcins_state->x[19];
	out_bias[2] = qcins_state->x[20];
	return true;
}

bool qcins_get_thrust(uintptr_t qcins_handle, float thrust[1])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	thrust[0] = qcins_state->x[21];
	return true;
}

bool qcins_get_mu(uintptr_t qcins_handle, float mu[1])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	mu[0] = qcins_state->x[22];
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
   const float x0 = x[4]*x[6];
   const float x1 = x[5]*x[7];
   const float x2 = expf(x[21])*param[0]*x[14];
   const float x3 = 2*x2;
   const float x4 = x[4]*x[7];
   const float x5 = x[5]*x[6];
   const float x6 = 2*x4 - 2*x5;
   const float x7 = expf(x[22]);
   const float x8 = x[4]*x[5];
   const float x9 = x[6]*x[7];
   const float x10 = x8 + x9;
   const float x11 = 2*x[3];
   const float x12 = (x[6]*x[6]);
   const float x13 = (x[5]*x[5]);
   const float x14 = -x13;
   const float x15 = (x[4]*x[4]);
   const float x16 = (x[7]*x[7]);
   const float x17 = x15 - x16;
   const float x18 = x12 + x14 + x17;
   const float x19 = x7*(x10*x11 + x18*x[2] - x6*x[1]);
   const float x20 = -x12;
   const float x21 = x13 + x17 + x20;
   const float x22 = 2*x4 + 2*x5;
   const float x23 = x0 - x1;
   const float x24 = x7*(-x11*x23 + x21*x[1] + x22*x[2]);
   const float x25 = (1.0F/2.0F)*x[8];
   const float x26 = (1.0F/2.0F)*x[6];
   const float x27 = (1.0F/2.0F)*x[7];
   const float x28 = (1.0F/2.0F)*x[4];
   const float x29 = (1.0F/2.0F)*x[5];
   const float x30 = Ts/param[5];
   xnew[0] = Ts*x[3] + x[0];
   xnew[1] = Ts*(x19*x6 - x21*x24 - x3*(x0 + x1)) + x[1];
   xnew[2] = Ts*(-x18*x19 - x22*x24 + x3*(x8 - x9)) + x[2];
   xnew[3] = Ts*(-2*x10*x19 - x2*(x14 + x15 + x16 + x20) + 2*x23*x24 + param[0]) + x[3];
   xnew[4] = Ts*(-x25*x[5] - x26*x[9] - x27*x[10]) + x[4];
   xnew[5] = Ts*(x25*x[4] + x26*x[10] - x27*x[9]) + x[5];
   xnew[6] = Ts*(x25*x[7] + x28*x[9] - x29*x[10]) + x[6];
   xnew[7] = Ts*(-x25*x[6] + x28*x[10] + x29*x[9]) + x[7];
   xnew[8] = Ts*param[1]*x[11] + x[8];
   xnew[9] = Ts*param[2]*x[12] + x[9];
   xnew[10] = Ts*(-(-u[2] + x[13] + x[17])*param[4] + param[3]*x[13]) + x[10];
   xnew[11] = x30*(u[0] - x[11] - x[15]) + x[11];
   xnew[12] = x30*(u[1] - x[12] - x[16]) + x[12];
   xnew[13] = x30*(u[2] - x[13] - x[17]) + x[13];
   xnew[14] = x30*(u[3] - x[14]) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = x[18];
   xnew[19] = x[19];
   xnew[20] = x[20];
   xnew[21] = x[21];
   xnew[22] = x[22];
}

void linearize_FH(const float * restrict x, const float Ts, const float * restrict param, float * restrict F, float * restrict H)
{
   const float x0 = expf(x[22]);
   const float x1 = 2*x[4];
   const float x2 = x1*x[7];
   const float x3 = 2*x[5];
   const float x4 = x3*x[6];
   const float x5 = -x2 + x4;
   const float x6 = (x[5]*x[5]);
   const float x7 = (x[6]*x[6]);
   const float x8 = -x7;
   const float x9 = (x[4]*x[4]);
   const float x10 = (x[7]*x[7]);
   const float x11 = -x10 + x9;
   const float x12 = x11 + x6 + x8;
   const float x13 = -x6;
   const float x14 = x11 + x13 + x7;
   const float x15 = x0*x5;
   const float x16 = x2 + x4;
   const float x17 = x0*x12;
   const float x18 = Ts*(-x14*x15 - x16*x17);
   const float x19 = x1*x[5];
   const float x20 = 2*x[7];
   const float x21 = x20*x[6];
   const float x22 = x19 + x21;
   const float x23 = x1*x[6];
   const float x24 = x3*x[7];
   const float x25 = -x23 + x24;
   const float x26 = Ts*(-x15*x22 - x17*x25);
   const float x27 = x1*x[2] - x20*x[1] + x3*x[3];
   const float x28 = x0*x27;
   const float x29 = x1*x[1];
   const float x30 = x20*x[2];
   const float x31 = 2*x[6];
   const float x32 = x31*x[3];
   const float x33 = x29 + x30 - x32;
   const float x34 = x12*x[1] + x16*x[2] + x25*x[3];
   const float x35 = x0*x34;
   const float x36 = x1*x35;
   const float x37 = x14*x[2] + x22*x[3] + x5*x[1];
   const float x38 = x0*x37;
   const float x39 = x20*x38;
   const float x40 = expf(x[21])*param[0];
   const float x41 = x40*x[14];
   const float x42 = x31*x41;
   const float x43 = -x36 + x39 - x42;
   const float x44 = x31*x[1];
   const float x45 = x3*x[2];
   const float x46 = x1*x[3];
   const float x47 = x44 - x45 + x46;
   const float x48 = x20*x[3] + x3*x[1] + x31*x[2];
   const float x49 = -x20*x41 - x3*x35 - x31*x38;
   const float x50 = x0*x48;
   const float x51 = -x44 + x45 - x46;
   const float x52 = x3*x38;
   const float x53 = x31*x35;
   const float x54 = x1*x41;
   const float x55 = -x52 + x53 - x54;
   const float x56 = x3*x41;
   const float x57 = -x29 - x30 + x32;
   const float x58 = x1*x38;
   const float x59 = x20*x35;
   const float x60 = x23 + x24;
   const float x61 = Ts*x40;
   const float x62 = Ts*x41;
   const float x63 = x0*x16;
   const float x64 = x0*x14;
   const float x65 = Ts*(-x22*x64 - x25*x63);
   const float x66 = x56 - x58 - x59;
   const float x67 = -x19 + x21;
   const float x68 = x0*x25;
   const float x69 = x0*x22;
   const float x70 = x10 + x13 + x8 + x9;
   const float x71 = (1.0F/2.0F)*Ts;
   const float x72 = x71*x[8];
   const float x73 = -x72;
   const float x74 = x71*x[9];
   const float x75 = -x74;
   const float x76 = x71*x[10];
   const float x77 = -x76;
   const float x78 = x71*x[5];
   const float x79 = -x78;
   const float x80 = x71*x[6];
   const float x81 = -x80;
   const float x82 = x71*x[7];
   const float x83 = -x82;
   const float x84 = x71*x[4];
   const float x85 = -Ts/param[5];
   const float x86 = -x50;
   const float x87 = -x28;
   const float x88 = x1*param[6];
   const float x89 = x20*param[7];
   const float x90 = x31*param[8];
   const float x91 = x88 + x89 - x90;
   const float x92 = x20*param[8] + x3*param[6] + x31*param[7];
   const float x93 = x31*param[6];
   const float x94 = x3*param[7];
   const float x95 = x1*param[8];
   const float x96 = x20*param[6];
   const float x97 = x1*param[7];
   const float x98 = x3*param[8];
   const float x99 = -x96 + x97 + x98;
   const float x100 = x93 - x94 + x95;
   F[0] = Ts;
   F[1] = Ts*(-x0*(x12*x12) - x0*(x5*x5));
   F[2] = x18;
   F[3] = x26;
   F[4] = Ts*(-x17*x33 - x28*x5 + x43);
   F[5] = Ts*(-x15*x47 - x17*x48 + x49);
   F[6] = Ts*(-x17*x51 - x5*x50 + x55);
   F[7] = Ts*(-x15*x57 - x17*x27 - x56 + x58 + x59);
   F[8] = -x60*x61;
   F[9] = -x60*x62;
   F[10] = Ts*(-x15*x37 - x17*x34);
   F[11] = x18;
   F[12] = Ts*(-x0*(x14*x14) - x0*(x16*x16));
   F[13] = x65;
   F[14] = Ts*(-x14*x28 - x33*x63 + x66);
   F[15] = Ts*(-x47*x64 - x48*x63 + x52 - x53 + x54);
   F[16] = Ts*(-x14*x50 + x49 - x51*x63);
   F[17] = Ts*(-x27*x63 + x43 - x57*x64);
   F[18] = -x61*x67;
   F[19] = -x62*x67;
   F[20] = Ts*(-x34*x63 - x37*x64);
   F[21] = x26;
   F[22] = x65;
   F[23] = Ts*(-x0*(x22*x22) - x0*(x25*x25));
   F[24] = Ts*(-x22*x28 - x33*x68 + x55);
   F[25] = Ts*(-x47*x69 - x48*x68 + x66);
   F[26] = Ts*(-x22*x50 + x36 - x39 + x42 - x51*x68);
   F[27] = Ts*(-x27*x68 + x49 - x57*x69);
   F[28] = -x61*x70;
   F[29] = -x62*x70;
   F[30] = Ts*(-x34*x68 - x37*x69);
   F[31] = x73;
   F[32] = x75;
   F[33] = x77;
   F[34] = x79;
   F[35] = x81;
   F[36] = x83;
   F[37] = x72;
   F[38] = x76;
   F[39] = x75;
   F[40] = x84;
   F[41] = x83;
   F[42] = x80;
   F[43] = x74;
   F[44] = x77;
   F[45] = x72;
   F[46] = x82;
   F[47] = x84;
   F[48] = x79;
   F[49] = x76;
   F[50] = x74;
   F[51] = x73;
   F[52] = x81;
   F[53] = x78;
   F[54] = x84;
   F[55] = Ts*param[1];
   F[56] = Ts*param[2];
   F[57] = Ts*(param[3] - param[4]);
   F[58] = -Ts*param[4];
   F[59] = x85;
   F[60] = x85;
   F[61] = x85;
   F[62] = x85;
   F[63] = x85;
   F[64] = x85;
   F[65] = x85;
   H[0] = 1;
   H[1] = -x17;
   H[2] = -x63;
   H[3] = -x68;
   H[4] = -x0*x33;
   H[5] = x86;
   H[6] = -x0*x51;
   H[7] = x87;
   H[8] = -x35;
   H[9] = -x15;
   H[10] = -x64;
   H[11] = -x69;
   H[12] = x87;
   H[13] = -x0*x47;
   H[14] = x86;
   H[15] = -x0*x57;
   H[16] = -x38;
   H[17] = -x40;
   H[18] = -x41;
   H[19] = 1;
   H[20] = 1;
   H[21] = 1;
   H[22] = 1;
   H[23] = 1;
   H[24] = 1;
   H[25] = x91;
   H[26] = x92;
   H[27] = -x93 + x94 - x95;
   H[28] = x99;
   H[29] = x99;
   H[30] = x100;
   H[31] = x92;
   H[32] = -x88 - x89 + x90;
   H[33] = x100;
   H[34] = x96 - x97 - x98;
   H[35] = x91;
   H[36] = x92;
}

void covariance_prediction(const float *restrict P, const float * restrict F, const float * restrict Q, float * restrict Pnew)
{
   const float x0 = F[0]*P[39] + P[3];
   const float x1 = F[0]*P[28] + P[2];
   const float x2 = F[0]*P[40] + P[4];
   const float x3 = F[0]*P[41] + P[5];
   const float x4 = F[0]*P[42] + P[6];
   const float x5 = F[0]*P[43] + P[7];
   const float x6 = F[0]*P[47] + P[11];
   const float x7 = F[0]*P[48] + P[12];
   const float x8 = F[0]*P[49] + P[13];
   const float x9 = F[1] + 1;
   const float x10 = F[0]*P[16] + P[1];
   const float x11 = F[12] + 1;
   const float x12 = F[23] + 1;
   const float x13 = F[0]*P[44] + P[8];
   const float x14 = F[0]*P[45] + P[9];
   const float x15 = F[0]*P[46] + P[10];
   const float x16 = F[65] + 1;
   const float x17 = x9*P[15] + F[2]*P[27] + F[3]*P[28] + F[4]*P[29] + F[5]*P[30] + F[6]*P[31] + F[7]*P[32] + F[8]*P[36] + F[9]*P[37] + F[10]*P[38];
   const float x18 = x9*P[16] + F[2]*P[28] + F[3]*P[39] + F[4]*P[40] + F[5]*P[41] + F[6]*P[42] + F[7]*P[43] + F[8]*P[47] + F[9]*P[48] + F[10]*P[49];
   const float x19 = x9*P[17] + F[2]*P[29] + F[3]*P[40] + F[4]*P[50] + F[5]*P[51] + F[6]*P[52] + F[7]*P[53] + F[8]*P[57] + F[9]*P[61] + F[10]*P[62];
   const float x20 = x9*P[18] + F[2]*P[30] + F[3]*P[41] + F[4]*P[51] + F[5]*P[63] + F[6]*P[64] + F[7]*P[65] + F[8]*P[69] + F[9]*P[73] + F[10]*P[74];
   const float x21 = x9*P[19] + F[2]*P[31] + F[3]*P[42] + F[4]*P[52] + F[5]*P[64] + F[6]*P[75] + F[7]*P[76] + F[8]*P[80] + F[9]*P[84] + F[10]*P[85];
   const float x22 = x9*P[20] + F[2]*P[32] + F[3]*P[43] + F[4]*P[53] + F[5]*P[65] + F[6]*P[76] + F[7]*P[86] + F[8]*P[90] + F[9]*P[94] + F[10]*P[95];
   const float x23 = x9*P[24] + F[2]*P[36] + F[3]*P[47] + F[4]*P[57] + F[5]*P[69] + F[6]*P[80] + F[7]*P[90] + F[8]*P[123] + F[9]*P[124] + F[10]*P[125];
   const float x24 = x9*P[25] + F[2]*P[37] + F[3]*P[48] + F[4]*P[61] + F[5]*P[73] + F[6]*P[84] + F[7]*P[94] + F[8]*P[124] + F[9]*P[141] + F[10]*P[142];
   const float x25 = x9*P[26] + F[2]*P[38] + F[3]*P[49] + F[4]*P[62] + F[5]*P[74] + F[6]*P[85] + F[7]*P[95] + F[8]*P[125] + F[9]*P[142] + F[10]*P[143];
   const float x26 = x9*P[14] + F[2]*P[15] + F[3]*P[16] + F[4]*P[17] + F[5]*P[18] + F[6]*P[19] + F[7]*P[20] + F[8]*P[24] + F[9]*P[25] + F[10]*P[26];
   const float x27 = x9*P[21] + F[2]*P[33] + F[3]*P[44] + F[4]*P[54] + F[5]*P[66] + F[6]*P[77] + F[7]*P[87] + F[10]*P[101];
   const float x28 = x9*P[22] + F[2]*P[34] + F[3]*P[45] + F[4]*P[55] + F[5]*P[67] + F[6]*P[78] + F[7]*P[88] + F[10]*P[105];
   const float x29 = x9*P[23] + F[2]*P[35] + F[3]*P[46] + F[4]*P[56] + F[5]*P[68] + F[6]*P[79] + F[7]*P[89] + F[10]*P[110];
   const float x30 = F[55]*P[114];
   const float x31 = F[56]*P[118];
   const float x32 = F[57]*P[122];
   const float x33 = F[58]*P[134];
   const float x34 = x11*P[15] + F[11]*P[14] + F[13]*P[16] + F[14]*P[17] + F[15]*P[18] + F[16]*P[19] + F[17]*P[20] + F[18]*P[24] + F[19]*P[25] + F[20]*P[26];
   const float x35 = x11*P[28] + F[11]*P[16] + F[13]*P[39] + F[14]*P[40] + F[15]*P[41] + F[16]*P[42] + F[17]*P[43] + F[18]*P[47] + F[19]*P[48] + F[20]*P[49];
   const float x36 = x11*P[29] + F[11]*P[17] + F[13]*P[40] + F[14]*P[50] + F[15]*P[51] + F[16]*P[52] + F[17]*P[53] + F[18]*P[57] + F[19]*P[61] + F[20]*P[62];
   const float x37 = x11*P[30] + F[11]*P[18] + F[13]*P[41] + F[14]*P[51] + F[15]*P[63] + F[16]*P[64] + F[17]*P[65] + F[18]*P[69] + F[19]*P[73] + F[20]*P[74];
   const float x38 = x11*P[31] + F[11]*P[19] + F[13]*P[42] + F[14]*P[52] + F[15]*P[64] + F[16]*P[75] + F[17]*P[76] + F[18]*P[80] + F[19]*P[84] + F[20]*P[85];
   const float x39 = x11*P[32] + F[11]*P[20] + F[13]*P[43] + F[14]*P[53] + F[15]*P[65] + F[16]*P[76] + F[17]*P[86] + F[18]*P[90] + F[19]*P[94] + F[20]*P[95];
   const float x40 = x11*P[36] + F[11]*P[24] + F[13]*P[47] + F[14]*P[57] + F[15]*P[69] + F[16]*P[80] + F[17]*P[90] + F[18]*P[123] + F[19]*P[124] + F[20]*P[125];
   const float x41 = x11*P[37] + F[11]*P[25] + F[13]*P[48] + F[14]*P[61] + F[15]*P[73] + F[16]*P[84] + F[17]*P[94] + F[18]*P[124] + F[19]*P[141] + F[20]*P[142];
   const float x42 = x11*P[38] + F[11]*P[26] + F[13]*P[49] + F[14]*P[62] + F[15]*P[74] + F[16]*P[85] + F[17]*P[95] + F[18]*P[125] + F[19]*P[142] + F[20]*P[143];
   const float x43 = x11*P[27] + F[11]*P[15] + F[13]*P[28] + F[14]*P[29] + F[15]*P[30] + F[16]*P[31] + F[17]*P[32] + F[18]*P[36] + F[19]*P[37] + F[20]*P[38];
   const float x44 = x11*P[33] + F[11]*P[21] + F[13]*P[44] + F[14]*P[54] + F[15]*P[66] + F[16]*P[77] + F[17]*P[87] + F[20]*P[101];
   const float x45 = x11*P[34] + F[11]*P[22] + F[13]*P[45] + F[14]*P[55] + F[15]*P[67] + F[16]*P[78] + F[17]*P[88] + F[20]*P[105];
   const float x46 = x11*P[35] + F[11]*P[23] + F[13]*P[46] + F[14]*P[56] + F[15]*P[68] + F[16]*P[79] + F[17]*P[89] + F[20]*P[110];
   const float x47 = x12*P[40] + F[21]*P[17] + F[22]*P[29] + F[24]*P[50] + F[25]*P[51] + F[26]*P[52] + F[27]*P[53] + F[28]*P[57] + F[29]*P[61] + F[30]*P[62];
   const float x48 = x12*P[41] + F[21]*P[18] + F[22]*P[30] + F[24]*P[51] + F[25]*P[63] + F[26]*P[64] + F[27]*P[65] + F[28]*P[69] + F[29]*P[73] + F[30]*P[74];
   const float x49 = x12*P[42] + F[21]*P[19] + F[22]*P[31] + F[24]*P[52] + F[25]*P[64] + F[26]*P[75] + F[27]*P[76] + F[28]*P[80] + F[29]*P[84] + F[30]*P[85];
   const float x50 = x12*P[43] + F[21]*P[20] + F[22]*P[32] + F[24]*P[53] + F[25]*P[65] + F[26]*P[76] + F[27]*P[86] + F[28]*P[90] + F[29]*P[94] + F[30]*P[95];
   const float x51 = x12*P[47] + F[21]*P[24] + F[22]*P[36] + F[24]*P[57] + F[25]*P[69] + F[26]*P[80] + F[27]*P[90] + F[28]*P[123] + F[29]*P[124] + F[30]*P[125];
   const float x52 = x12*P[48] + F[21]*P[25] + F[22]*P[37] + F[24]*P[61] + F[25]*P[73] + F[26]*P[84] + F[27]*P[94] + F[28]*P[124] + F[29]*P[141] + F[30]*P[142];
   const float x53 = x12*P[49] + F[21]*P[26] + F[22]*P[38] + F[24]*P[62] + F[25]*P[74] + F[26]*P[85] + F[27]*P[95] + F[28]*P[125] + F[29]*P[142] + F[30]*P[143];
   const float x54 = x12*P[44] + F[21]*P[21] + F[22]*P[33] + F[24]*P[54] + F[25]*P[66] + F[26]*P[77] + F[27]*P[87] + F[30]*P[101];
   const float x55 = x12*P[45] + F[21]*P[22] + F[22]*P[34] + F[24]*P[55] + F[25]*P[67] + F[26]*P[78] + F[27]*P[88] + F[30]*P[105];
   const float x56 = x12*P[46] + F[21]*P[23] + F[22]*P[35] + F[24]*P[56] + F[25]*P[68] + F[26]*P[79] + F[27]*P[89] + F[30]*P[110];
   const float x57 = F[31]*P[66] + F[32]*P[77] + F[33]*P[87] + F[34]*P[96] + P[54];
   const float x58 = F[31]*P[67] + F[32]*P[78] + F[33]*P[88] + F[35]*P[102] + P[55];
   const float x59 = F[31]*P[68] + F[32]*P[79] + F[33]*P[89] + F[36]*P[106] + P[56];
   const float x60 = F[31]*P[63] + F[32]*P[64] + F[33]*P[65] + F[34]*P[66] + F[35]*P[67] + F[36]*P[68] + P[51];
   const float x61 = F[31]*P[64] + F[32]*P[75] + F[33]*P[76] + F[34]*P[77] + F[35]*P[78] + F[36]*P[79] + P[52];
   const float x62 = F[31]*P[65] + F[32]*P[76] + F[33]*P[86] + F[34]*P[87] + F[35]*P[88] + F[36]*P[89] + P[53];
   const float x63 = F[31]*P[51] + F[32]*P[52] + F[33]*P[53] + F[34]*P[54] + F[35]*P[55] + F[36]*P[56] + P[50];
   const float x64 = F[55]*P[97];
   const float x65 = F[56]*P[98];
   const float x66 = F[57]*P[107];
   const float x67 = F[58]*P[108];
   const float x68 = F[37]*P[54] + F[38]*P[77] + F[39]*P[87] + F[40]*P[96] + P[66];
   const float x69 = F[37]*P[55] + F[38]*P[78] + F[39]*P[88] + F[41]*P[102] + P[67];
   const float x70 = F[37]*P[56] + F[38]*P[79] + F[39]*P[89] + F[42]*P[106] + P[68];
   const float x71 = F[37]*P[50] + F[38]*P[52] + F[39]*P[53] + F[40]*P[54] + F[41]*P[55] + F[42]*P[56] + P[51];
   const float x72 = F[37]*P[52] + F[38]*P[75] + F[39]*P[76] + F[40]*P[77] + F[41]*P[78] + F[42]*P[79] + P[64];
   const float x73 = F[37]*P[53] + F[38]*P[76] + F[39]*P[86] + F[40]*P[87] + F[41]*P[88] + F[42]*P[89] + P[65];
   const float x74 = F[37]*P[51] + F[38]*P[64] + F[39]*P[65] + F[40]*P[66] + F[41]*P[67] + F[42]*P[68] + P[63];
   const float x75 = F[43]*P[54] + F[44]*P[66] + F[45]*P[87] + F[46]*P[96] + P[77];
   const float x76 = F[43]*P[55] + F[44]*P[67] + F[45]*P[88] + F[47]*P[102] + P[78];
   const float x77 = F[43]*P[56] + F[44]*P[68] + F[45]*P[89] + F[48]*P[106] + P[79];
   const float x78 = F[43]*P[50] + F[44]*P[51] + F[45]*P[53] + F[46]*P[54] + F[47]*P[55] + F[48]*P[56] + P[52];
   const float x79 = F[43]*P[51] + F[44]*P[63] + F[45]*P[65] + F[46]*P[66] + F[47]*P[67] + F[48]*P[68] + P[64];
   const float x80 = F[43]*P[53] + F[44]*P[65] + F[45]*P[86] + F[46]*P[87] + F[47]*P[88] + F[48]*P[89] + P[76];
   const float x81 = F[43]*P[52] + F[44]*P[64] + F[45]*P[76] + F[46]*P[77] + F[47]*P[78] + F[48]*P[79] + P[75];
   const float x82 = F[49]*P[54] + F[50]*P[66] + F[51]*P[77] + F[52]*P[96] + P[87];
   const float x83 = F[49]*P[55] + F[50]*P[67] + F[51]*P[78] + F[53]*P[102] + P[88];
   const float x84 = F[49]*P[56] + F[50]*P[68] + F[51]*P[79] + F[54]*P[106] + P[89];
   const float x85 = F[55]*P[111] + P[97];
   const float x86 = F[55]*P[112] + P[99];
   const float x87 = F[59] + 1;
   const float x88 = F[61] + 1;
   const float x89 = F[57]*P[119] + F[58]*P[120] + P[107];
   const float x90 = F[57]*P[120] + F[58]*P[132] + P[108];
   const float x91 = F[63] + 1;
   const float x92 = x87*P[112] + F[60]*P[126];
   const float x93 = x88*P[116] + F[62]*P[129];
   const float x94 = x91*P[120] + F[64]*P[132];
   Pnew[0] = x0*F[0] + F[0]*P[3] + P[0] + Q[0];
   Pnew[1] = x0*F[3] + x1*F[2] + x10*x9 + x2*F[4] + x3*F[5] + x4*F[6] + x5*F[7] + x6*F[8] + x7*F[9] + x8*F[10];
   Pnew[2] = x0*F[13] + x1*x11 + x10*F[11] + x2*F[14] + x3*F[15] + x4*F[16] + x5*F[17] + x6*F[18] + x7*F[19] + x8*F[20];
   Pnew[3] = x0*x12 + x1*F[22] + x10*F[21] + x2*F[24] + x3*F[25] + x4*F[26] + x5*F[27] + x6*F[28] + x7*F[29] + x8*F[30];
   Pnew[4] = x13*F[34] + x14*F[35] + x15*F[36] + x2 + x3*F[31] + x4*F[32] + x5*F[33];
   Pnew[5] = x13*F[40] + x14*F[41] + x15*F[42] + x2*F[37] + x3 + x4*F[38] + x5*F[39];
   Pnew[6] = x13*F[46] + x14*F[47] + x15*F[48] + x2*F[43] + x3*F[44] + x4 + x5*F[45];
   Pnew[7] = x13*F[52] + x14*F[53] + x15*F[54] + x2*F[49] + x3*F[50] + x4*F[51] + x5;
   Pnew[8] = x13;
   Pnew[9] = x14;
   Pnew[10] = x15;
   Pnew[11] = x16*x6;
   Pnew[12] = x7;
   Pnew[13] = x8;
   Pnew[14] = x17*F[2] + x18*F[3] + x19*F[4] + x20*F[5] + x21*F[6] + x22*F[7] + x23*F[8] + x24*F[9] + x25*F[10] + x26*x9 + Q[1];
   Pnew[15] = x11*x17 + x18*F[13] + x19*F[14] + x20*F[15] + x21*F[16] + x22*F[17] + x23*F[18] + x24*F[19] + x25*F[20] + x26*F[11];
   Pnew[16] = x12*x18 + x17*F[22] + x19*F[24] + x20*F[25] + x21*F[26] + x22*F[27] + x23*F[28] + x24*F[29] + x25*F[30] + x26*F[21];
   Pnew[17] = x19 + x20*F[31] + x21*F[32] + x22*F[33] + x27*F[34] + x28*F[35] + x29*F[36];
   Pnew[18] = x19*F[37] + x20 + x21*F[38] + x22*F[39] + x27*F[40] + x28*F[41] + x29*F[42];
   Pnew[19] = x19*F[43] + x20*F[44] + x21 + x22*F[45] + x27*F[46] + x28*F[47] + x29*F[48];
   Pnew[20] = x19*F[49] + x20*F[50] + x21*F[51] + x22 + x27*F[52] + x28*F[53] + x29*F[54];
   Pnew[21] = x27 + x30*F[10];
   Pnew[22] = x28 + x31*F[10];
   Pnew[23] = x29 + x32*F[10] + x33*F[10];
   Pnew[24] = x16*x23;
   Pnew[25] = x24;
   Pnew[26] = x25;
   Pnew[27] = x11*x43 + x34*F[11] + x35*F[13] + x36*F[14] + x37*F[15] + x38*F[16] + x39*F[17] + x40*F[18] + x41*F[19] + x42*F[20] + Q[2];
   Pnew[28] = x12*x35 + x34*F[21] + x36*F[24] + x37*F[25] + x38*F[26] + x39*F[27] + x40*F[28] + x41*F[29] + x42*F[30] + x43*F[22];
   Pnew[29] = x36 + x37*F[31] + x38*F[32] + x39*F[33] + x44*F[34] + x45*F[35] + x46*F[36];
   Pnew[30] = x36*F[37] + x37 + x38*F[38] + x39*F[39] + x44*F[40] + x45*F[41] + x46*F[42];
   Pnew[31] = x36*F[43] + x37*F[44] + x38 + x39*F[45] + x44*F[46] + x45*F[47] + x46*F[48];
   Pnew[32] = x36*F[49] + x37*F[50] + x38*F[51] + x39 + x44*F[52] + x45*F[53] + x46*F[54];
   Pnew[33] = x30*F[20] + x44;
   Pnew[34] = x31*F[20] + x45;
   Pnew[35] = x32*F[20] + x33*F[20] + x46;
   Pnew[36] = x16*x40;
   Pnew[37] = x41;
   Pnew[38] = x42;
   Pnew[39] = x12*(x12*P[39] + F[21]*P[16] + F[22]*P[28] + F[24]*P[40] + F[25]*P[41] + F[26]*P[42] + F[27]*P[43] + F[28]*P[47] + F[29]*P[48] + F[30]*P[49]) + x47*F[24] + x48*F[25] + x49*F[26] + x50*F[27] + x51*F[28] + x52*F[29] + x53*F[30] + (x12*P[16] + F[21]*P[14] + F[22]*P[15] + F[24]*P[17] + F[25]*P[18] + F[26]*P[19] + F[27]*P[20] + F[28]*P[24] + F[29]*P[25] + F[30]*P[26])*F[21] + (x12*P[28] + F[21]*P[15] + F[22]*P[27] + F[24]*P[29] + F[25]*P[30] + F[26]*P[31] + F[27]*P[32] + F[28]*P[36] + F[29]*P[37] + F[30]*P[38])*F[22] + Q[3];
   Pnew[40] = x47 + x48*F[31] + x49*F[32] + x50*F[33] + x54*F[34] + x55*F[35] + x56*F[36];
   Pnew[41] = x47*F[37] + x48 + x49*F[38] + x50*F[39] + x54*F[40] + x55*F[41] + x56*F[42];
   Pnew[42] = x47*F[43] + x48*F[44] + x49 + x50*F[45] + x54*F[46] + x55*F[47] + x56*F[48];
   Pnew[43] = x47*F[49] + x48*F[50] + x49*F[51] + x50 + x54*F[52] + x55*F[53] + x56*F[54];
   Pnew[44] = x30*F[30] + x54;
   Pnew[45] = x31*F[30] + x55;
   Pnew[46] = x32*F[30] + x33*F[30] + x56;
   Pnew[47] = x16*x51;
   Pnew[48] = x52;
   Pnew[49] = x53;
   Pnew[50] = x57*F[34] + x58*F[35] + x59*F[36] + x60*F[31] + x61*F[32] + x62*F[33] + x63 + Q[4];
   Pnew[51] = x57*F[40] + x58*F[41] + x59*F[42] + x60 + x61*F[38] + x62*F[39] + x63*F[37];
   Pnew[52] = x57*F[46] + x58*F[47] + x59*F[48] + x60*F[44] + x61 + x62*F[45] + x63*F[43];
   Pnew[53] = x57*F[52] + x58*F[53] + x59*F[54] + x60*F[50] + x61*F[51] + x62 + x63*F[49];
   Pnew[54] = x57 + x64*F[34];
   Pnew[55] = x58 + x65*F[34];
   Pnew[56] = x59 + x66*F[36] + x67*F[36];
   Pnew[57] = x16*(F[31]*P[69] + F[32]*P[80] + F[33]*P[90] + P[57]);
   Pnew[58] = F[31]*P[70] + F[32]*P[81] + F[33]*P[91] + F[34]*P[100] + P[58];
   Pnew[59] = F[31]*P[71] + F[32]*P[82] + F[33]*P[92] + F[35]*P[104] + P[59];
   Pnew[60] = F[31]*P[72] + F[32]*P[83] + F[33]*P[93] + F[36]*P[109] + P[60];
   Pnew[61] = F[31]*P[73] + F[32]*P[84] + F[33]*P[94] + P[61];
   Pnew[62] = F[31]*P[74] + F[32]*P[85] + F[33]*P[95] + F[34]*P[101] + F[35]*P[105] + F[36]*P[110] + P[62];
   Pnew[63] = x68*F[40] + x69*F[41] + x70*F[42] + x71*F[37] + x72*F[38] + x73*F[39] + x74 + Q[5];
   Pnew[64] = x68*F[46] + x69*F[47] + x70*F[48] + x71*F[43] + x72 + x73*F[45] + x74*F[44];
   Pnew[65] = x68*F[52] + x69*F[53] + x70*F[54] + x71*F[49] + x72*F[51] + x73 + x74*F[50];
   Pnew[66] = x64*F[40] + x68;
   Pnew[67] = x65*F[40] + x69;
   Pnew[68] = x66*F[42] + x67*F[42] + x70;
   Pnew[69] = x16*(F[37]*P[57] + F[38]*P[80] + F[39]*P[90] + P[69]);
   Pnew[70] = F[37]*P[58] + F[38]*P[81] + F[39]*P[91] + F[40]*P[100] + P[70];
   Pnew[71] = F[37]*P[59] + F[38]*P[82] + F[39]*P[92] + F[41]*P[104] + P[71];
   Pnew[72] = F[37]*P[60] + F[38]*P[83] + F[39]*P[93] + F[42]*P[109] + P[72];
   Pnew[73] = F[37]*P[61] + F[38]*P[84] + F[39]*P[94] + P[73];
   Pnew[74] = F[37]*P[62] + F[38]*P[85] + F[39]*P[95] + F[40]*P[101] + F[41]*P[105] + F[42]*P[110] + P[74];
   Pnew[75] = x75*F[46] + x76*F[47] + x77*F[48] + x78*F[43] + x79*F[44] + x80*F[45] + x81 + Q[6];
   Pnew[76] = x75*F[52] + x76*F[53] + x77*F[54] + x78*F[49] + x79*F[50] + x80 + x81*F[51];
   Pnew[77] = x64*F[46] + x75;
   Pnew[78] = x65*F[46] + x76;
   Pnew[79] = x66*F[48] + x67*F[48] + x77;
   Pnew[80] = x16*(F[43]*P[57] + F[44]*P[69] + F[45]*P[90] + P[80]);
   Pnew[81] = F[43]*P[58] + F[44]*P[70] + F[45]*P[91] + F[46]*P[100] + P[81];
   Pnew[82] = F[43]*P[59] + F[44]*P[71] + F[45]*P[92] + F[47]*P[104] + P[82];
   Pnew[83] = F[43]*P[60] + F[44]*P[72] + F[45]*P[93] + F[48]*P[109] + P[83];
   Pnew[84] = F[43]*P[61] + F[44]*P[73] + F[45]*P[94] + P[84];
   Pnew[85] = F[43]*P[62] + F[44]*P[74] + F[45]*P[95] + F[46]*P[101] + F[47]*P[105] + F[48]*P[110] + P[85];
   Pnew[86] = x82*F[52] + x83*F[53] + x84*F[54] + (F[49]*P[50] + F[50]*P[51] + F[51]*P[52] + F[52]*P[54] + F[53]*P[55] + F[54]*P[56] + P[53])*F[49] + (F[49]*P[51] + F[50]*P[63] + F[51]*P[64] + F[52]*P[66] + F[53]*P[67] + F[54]*P[68] + P[65])*F[50] + (F[49]*P[52] + F[50]*P[64] + F[51]*P[75] + F[52]*P[77] + F[53]*P[78] + F[54]*P[79] + P[76])*F[51] + F[49]*P[53] + F[50]*P[65] + F[51]*P[76] + F[52]*P[87] + F[53]*P[88] + F[54]*P[89] + P[86] + Q[7];
   Pnew[87] = x64*F[52] + x82;
   Pnew[88] = x65*F[52] + x83;
   Pnew[89] = x66*F[54] + x67*F[54] + x84;
   Pnew[90] = x16*(F[49]*P[57] + F[50]*P[69] + F[51]*P[80] + P[90]);
   Pnew[91] = F[49]*P[58] + F[50]*P[70] + F[51]*P[81] + F[52]*P[100] + P[91];
   Pnew[92] = F[49]*P[59] + F[50]*P[71] + F[51]*P[82] + F[53]*P[104] + P[92];
   Pnew[93] = F[49]*P[60] + F[50]*P[72] + F[51]*P[83] + F[54]*P[109] + P[93];
   Pnew[94] = F[49]*P[61] + F[50]*P[73] + F[51]*P[84] + P[94];
   Pnew[95] = F[49]*P[62] + F[50]*P[74] + F[51]*P[85] + F[52]*P[101] + F[53]*P[105] + F[54]*P[110] + P[95];
   Pnew[96] = x64 + x85*F[55] + P[96] + Q[8];
   Pnew[97] = x85*x87 + x86*F[60];
   Pnew[98] = x88*P[98];
   Pnew[99] = x86;
   Pnew[100] = F[55]*P[113] + P[100];
   Pnew[101] = x30 + P[101];
   Pnew[102] = (F[56]*F[56])*P[115] + P[102] + Q[9];
   Pnew[103] = F[56]*P[116] + P[103];
   Pnew[104] = F[56]*P[117] + P[104];
   Pnew[105] = x31 + P[105];
   Pnew[106] = x66 + x67 + x89*F[57] + x90*F[58] + P[106] + Q[10];
   Pnew[107] = x89*x91 + x90*F[64];
   Pnew[108] = x90;
   Pnew[109] = F[57]*P[121] + F[58]*P[133] + P[109];
   Pnew[110] = x32 + x33 + P[110];
   Pnew[111] = x87*(x87*P[111] + F[60]*P[112]) + x92*F[60] + Q[11];
   Pnew[112] = x92;
   Pnew[113] = x87*P[113] + F[60]*P[127];
   Pnew[114] = x87*P[114] + F[60]*P[128];
   Pnew[115] = x88*(x88*P[115] + F[62]*P[116]) + x93*F[62] + Q[12];
   Pnew[116] = x93;
   Pnew[117] = x88*P[117] + F[62]*P[130];
   Pnew[118] = x88*P[118] + F[62]*P[131];
   Pnew[119] = x91*(x91*P[119] + F[64]*P[120]) + x94*F[64] + Q[13];
   Pnew[120] = x94;
   Pnew[121] = x91*P[121] + F[64]*P[133];
   Pnew[122] = x91*P[122] + F[64]*P[134];
   Pnew[123] = (x16*x16)*P[123] + Q[14];
   Pnew[124] = x16*P[124];
   Pnew[125] = x16*P[125];
   Pnew[126] = P[126] + Q[15];
   Pnew[127] = P[127];
   Pnew[128] = P[128];
   Pnew[129] = P[129] + Q[16];
   Pnew[130] = P[130];
   Pnew[131] = P[131];
   Pnew[132] = P[132] + Q[17];
   Pnew[133] = P[133];
   Pnew[134] = P[134];
   Pnew[135] = P[135] + Q[18];
   Pnew[136] = P[136];
   Pnew[137] = P[137] + Q[19];
   Pnew[138] = P[138];
   Pnew[139] = P[139] + Q[20];
   Pnew[140] = P[140];
   Pnew[141] = P[141] + Q[21];
   Pnew[142] = P[142];
   Pnew[143] = P[143] + Q[22];
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
   const float x6 = x5*P[1];
   const float x7 = x5*P[2];
   const float x8 = x5*P[3];
   const float x9 = x5*P[4];
   const float x10 = x5*P[5];
   const float x11 = x5*P[6];
   const float x12 = x5*P[7];
   const float x13 = x5*P[13];
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
   xnew[18] = x[18];
   xnew[19] = x[19];
   xnew[20] = x[20];
   xnew[21] = -x3*P[12] + x[21];
   xnew[22] = -x3*P[13] + x[22];
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
   Pnew[13] = x4*P[13];
   Pnew[14] = -x5*(P[1]*P[1]) + P[14];
   Pnew[15] = -x6*P[2] + P[15];
   Pnew[16] = -x6*P[3] + P[16];
   Pnew[17] = -x6*P[4] + P[17];
   Pnew[18] = -x6*P[5] + P[18];
   Pnew[19] = -x6*P[6] + P[19];
   Pnew[20] = -x6*P[7] + P[20];
   Pnew[21] = -x6*P[8] + P[21];
   Pnew[22] = -x6*P[9] + P[22];
   Pnew[23] = -x6*P[10] + P[23];
   Pnew[24] = -x6*P[11] + P[24];
   Pnew[25] = -x6*P[12] + P[25];
   Pnew[26] = -x6*P[13] + P[26];
   Pnew[27] = -x5*(P[2]*P[2]) + P[27];
   Pnew[28] = -x7*P[3] + P[28];
   Pnew[29] = -x7*P[4] + P[29];
   Pnew[30] = -x7*P[5] + P[30];
   Pnew[31] = -x7*P[6] + P[31];
   Pnew[32] = -x7*P[7] + P[32];
   Pnew[33] = -x7*P[8] + P[33];
   Pnew[34] = -x7*P[9] + P[34];
   Pnew[35] = -x7*P[10] + P[35];
   Pnew[36] = -x7*P[11] + P[36];
   Pnew[37] = -x7*P[12] + P[37];
   Pnew[38] = -x7*P[13] + P[38];
   Pnew[39] = -x5*(P[3]*P[3]) + P[39];
   Pnew[40] = -x8*P[4] + P[40];
   Pnew[41] = -x8*P[5] + P[41];
   Pnew[42] = -x8*P[6] + P[42];
   Pnew[43] = -x8*P[7] + P[43];
   Pnew[44] = -x8*P[8] + P[44];
   Pnew[45] = -x8*P[9] + P[45];
   Pnew[46] = -x8*P[10] + P[46];
   Pnew[47] = -x8*P[11] + P[47];
   Pnew[48] = -x8*P[12] + P[48];
   Pnew[49] = -x8*P[13] + P[49];
   Pnew[50] = -x5*(P[4]*P[4]) + P[50];
   Pnew[51] = -x9*P[5] + P[51];
   Pnew[52] = -x9*P[6] + P[52];
   Pnew[53] = -x9*P[7] + P[53];
   Pnew[54] = -x9*P[8] + P[54];
   Pnew[55] = -x9*P[9] + P[55];
   Pnew[56] = -x9*P[10] + P[56];
   Pnew[57] = -x9*P[11] + P[57];
   Pnew[58] = P[58];
   Pnew[59] = P[59];
   Pnew[60] = P[60];
   Pnew[61] = -x9*P[12] + P[61];
   Pnew[62] = -x9*P[13] + P[62];
   Pnew[63] = -x5*(P[5]*P[5]) + P[63];
   Pnew[64] = -x10*P[6] + P[64];
   Pnew[65] = -x10*P[7] + P[65];
   Pnew[66] = -x10*P[8] + P[66];
   Pnew[67] = -x10*P[9] + P[67];
   Pnew[68] = -x10*P[10] + P[68];
   Pnew[69] = -x10*P[11] + P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = P[72];
   Pnew[73] = -x10*P[12] + P[73];
   Pnew[74] = -x10*P[13] + P[74];
   Pnew[75] = -x5*(P[6]*P[6]) + P[75];
   Pnew[76] = -x11*P[7] + P[76];
   Pnew[77] = -x11*P[8] + P[77];
   Pnew[78] = -x11*P[9] + P[78];
   Pnew[79] = -x11*P[10] + P[79];
   Pnew[80] = -x11*P[11] + P[80];
   Pnew[81] = P[81];
   Pnew[82] = P[82];
   Pnew[83] = P[83];
   Pnew[84] = -x11*P[12] + P[84];
   Pnew[85] = -x11*P[13] + P[85];
   Pnew[86] = -x5*(P[7]*P[7]) + P[86];
   Pnew[87] = -x12*P[8] + P[87];
   Pnew[88] = -x12*P[9] + P[88];
   Pnew[89] = -x12*P[10] + P[89];
   Pnew[90] = -x12*P[11] + P[90];
   Pnew[91] = P[91];
   Pnew[92] = P[92];
   Pnew[93] = P[93];
   Pnew[94] = -x12*P[12] + P[94];
   Pnew[95] = -x12*P[13] + P[95];
   Pnew[96] = -x5*(P[8]*P[8]) + P[96];
   Pnew[97] = P[97];
   Pnew[98] = P[98];
   Pnew[99] = P[99];
   Pnew[100] = P[100];
   Pnew[101] = -x13*P[8] + P[101];
   Pnew[102] = -x5*(P[9]*P[9]) + P[102];
   Pnew[103] = P[103];
   Pnew[104] = P[104];
   Pnew[105] = -x13*P[9] + P[105];
   Pnew[106] = -x5*(P[10]*P[10]) + P[106];
   Pnew[107] = P[107];
   Pnew[108] = P[108];
   Pnew[109] = P[109];
   Pnew[110] = -x13*P[10] + P[110];
   Pnew[111] = P[111];
   Pnew[112] = P[112];
   Pnew[113] = P[113];
   Pnew[114] = P[114];
   Pnew[115] = P[115];
   Pnew[116] = P[116];
   Pnew[117] = P[117];
   Pnew[118] = P[118];
   Pnew[119] = P[119];
   Pnew[120] = P[120];
   Pnew[121] = P[121];
   Pnew[122] = P[122];
   Pnew[123] = -x5*(P[11]*P[11]) + P[123];
   Pnew[124] = -x5*P[11]*P[12] + P[124];
   Pnew[125] = -x13*P[11] + P[125];
   Pnew[126] = P[126];
   Pnew[127] = P[127];
   Pnew[128] = P[128];
   Pnew[129] = P[129];
   Pnew[130] = P[130];
   Pnew[131] = P[131];
   Pnew[132] = P[132];
   Pnew[133] = P[133];
   Pnew[134] = P[134];
   Pnew[135] = P[135];
   Pnew[136] = P[136];
   Pnew[137] = P[137];
   Pnew[138] = P[138];
   Pnew[139] = P[139];
   Pnew[140] = P[140];
   Pnew[141] = -x5*(P[12]*P[12]) + P[141];
   Pnew[142] = -x13*P[12] + P[142];
   Pnew[143] = -x5*(P[13]*P[13]) + P[143];
}

void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, const float * restrict param,
        float *restrict xnew, float *restrict Pnew)
{
   const float x0 = 2*param[6];
   const float x1 = x[4]*x[6];
   const float x2 = x[5]*x[7];
   const float x3 = 2*param[7];
   const float x4 = x[4]*x[5];
   const float x5 = x[6]*x[7];
   const float x6 = (x[4]*x[4]);
   const float x7 = (x[5]*x[5]);
   const float x8 = x6 - x7;
   const float x9 = (x[7]*x[7]);
   const float x10 = (x[6]*x[6]);
   const float x11 = -x10;
   const float x12 = x0*(x1 + x2) - x3*(x4 - x5) + (x11 + x8 + x9)*param[8] - mag[2];
   const float x13 = H[25]*P[50] + H[26]*P[51] + H[27]*P[52] + H[28]*P[53];
   const float x14 = H[25]*P[51] + H[26]*P[63] + H[27]*P[64] + H[28]*P[65];
   const float x15 = H[25]*P[52] + H[26]*P[64] + H[27]*P[75] + H[28]*P[76];
   const float x16 = H[25]*P[53] + H[26]*P[65] + H[27]*P[76] + H[28]*P[86];
   const float x17 = x13*H[29] + x14*H[30] + x15*H[31] + x16*H[32];
   const float x18 = H[29]*P[50] + H[30]*P[51] + H[31]*P[52] + H[32]*P[53];
   const float x19 = H[29]*P[51] + H[30]*P[63] + H[31]*P[64] + H[32]*P[65];
   const float x20 = H[29]*P[52] + H[30]*P[64] + H[31]*P[75] + H[32]*P[76];
   const float x21 = H[29]*P[53] + H[30]*P[65] + H[31]*P[76] + H[32]*P[86];
   const float x22 = x18*H[33] + x19*H[34] + x20*H[35] + x21*H[36];
   const float x23 = x17*x22;
   const float x24 = x13*H[33] + x14*H[34] + x15*H[35] + x16*H[36];
   const float x25 = x18*H[29] + x19*H[30] + x20*H[31] + x21*H[32] + R[8];
   const float x26 = H[33]*P[50] + H[34]*P[51] + H[35]*P[52] + H[36]*P[53];
   const float x27 = H[33]*P[51] + H[34]*P[63] + H[35]*P[64] + H[36]*P[65];
   const float x28 = H[33]*P[52] + H[34]*P[64] + H[35]*P[75] + H[36]*P[76];
   const float x29 = H[33]*P[53] + H[34]*P[65] + H[35]*P[76] + H[36]*P[86];
   const float x30 = x26*H[25] + x27*H[26] + x28*H[27] + x29*H[28];
   const float x31 = x26*H[29] + x27*H[30] + x28*H[31] + x29*H[32];
   const float x32 = x18*H[25] + x19*H[26] + x20*H[27] + x21*H[28];
   const float x33 = x24*x32;
   const float x34 = x26*H[33] + x27*H[34] + x28*H[35] + x29*H[36] + R[9];
   const float x35 = x17*x32;
   const float x36 = x24*x30;
   const float x37 = x13*H[25] + x14*H[26] + x15*H[27] + x16*H[28] + R[7];
   const float x38 = x22*x31;
   const float x39 = x25*x37;
   const float x40 = 1.0F/(x23*x30 - x25*x36 + x31*x33 - x34*x35 + x34*x39 - x37*x38);
   const float x41 = x40*(x23 - x24*x25);
   const float x42 = x13*x41;
   const float x43 = x40*(-x22*x37 + x33);
   const float x44 = x18*x43;
   const float x45 = x40*(-x35 + x39);
   const float x46 = x26*x45;
   const float x47 = 2*param[8];
   const float x48 = x[4]*x[7];
   const float x49 = x[5]*x[6];
   const float x50 = -x9;
   const float x51 = -x0*(x48 - x49) + x47*(x4 + x5) + (x10 + x50 + x8)*param[7] - mag[1];
   const float x52 = x40*(x17*x30 - x31*x37);
   const float x53 = x26*x52;
   const float x54 = x40*(-x17*x34 + x24*x31);
   const float x55 = x13*x54;
   const float x56 = x40*(x34*x37 - x36);
   const float x57 = x18*x56;
   const float x58 = x3*(x48 + x49) - x47*(x1 - x2) + (x11 + x50 + x6 + x7)*param[6] - mag[0];
   const float x59 = x40*(-x25*x30 + x31*x32);
   const float x60 = x26*x59;
   const float x61 = x40*(x22*x30 - x32*x34);
   const float x62 = x18*x61;
   const float x63 = x40*(x25*x34 - x38);
   const float x64 = x13*x63;
   const float x65 = x14*x41;
   const float x66 = x19*x43;
   const float x67 = x27*x45;
   const float x68 = x27*x52;
   const float x69 = x14*x54;
   const float x70 = x19*x56;
   const float x71 = x27*x59;
   const float x72 = x19*x61;
   const float x73 = x14*x63;
   const float x74 = x15*x41;
   const float x75 = x20*x43;
   const float x76 = x28*x45;
   const float x77 = x28*x52;
   const float x78 = x15*x54;
   const float x79 = x20*x56;
   const float x80 = x28*x59;
   const float x81 = x20*x61;
   const float x82 = x15*x63;
   const float x83 = x16*x41;
   const float x84 = x21*x43;
   const float x85 = x29*x45;
   const float x86 = x29*x52;
   const float x87 = x16*x54;
   const float x88 = x21*x56;
   const float x89 = x29*x59;
   const float x90 = x21*x61;
   const float x91 = x16*x63;
   const float x92 = H[25]*P[56] + H[26]*P[68] + H[27]*P[79] + H[28]*P[89];
   const float x93 = x41*x92;
   const float x94 = H[29]*P[56] + H[30]*P[68] + H[31]*P[79] + H[32]*P[89];
   const float x95 = x43*x94;
   const float x96 = H[33]*P[56] + H[34]*P[68] + H[35]*P[79] + H[36]*P[89];
   const float x97 = x45*x96;
   const float x98 = x52*x96;
   const float x99 = x54*x92;
   const float x100 = x56*x94;
   const float x101 = x59*x96;
   const float x102 = x61*x94;
   const float x103 = x63*x92;
   const float x104 = H[25]*P[60] + H[26]*P[72] + H[27]*P[83] + H[28]*P[93];
   const float x105 = x104*x41;
   const float x106 = H[29]*P[60] + H[30]*P[72] + H[31]*P[83] + H[32]*P[93];
   const float x107 = x106*x43;
   const float x108 = H[33]*P[60] + H[34]*P[72] + H[35]*P[83] + H[36]*P[93];
   const float x109 = x108*x45;
   const float x110 = x108*x52;
   const float x111 = x104*x54;
   const float x112 = x106*x56;
   const float x113 = x108*x59;
   const float x114 = x106*x61;
   const float x115 = x104*x63;
   const float x116 = x42 + x44 + x46;
   const float x117 = x53 + x55 + x57;
   const float x118 = x60 + x62 + x64;
   const float x119 = -x116*H[34] - x117*H[30] - x118*H[26];
   const float x120 = -x116*H[35] - x117*H[31] - x118*H[27];
   const float x121 = -x116*H[36] - x117*H[32] - x118*H[28];
   const float x122 = -x116*H[33] - x117*H[29] - x118*H[25] + 1;
   const float x123 = x65 + x66 + x67;
   const float x124 = x68 + x69 + x70;
   const float x125 = x71 + x72 + x73;
   const float x126 = -x123*H[33] - x124*H[29] - x125*H[25];
   const float x127 = -x123*H[35] - x124*H[31] - x125*H[27];
   const float x128 = -x123*H[36] - x124*H[32] - x125*H[28];
   const float x129 = -x123*H[34] - x124*H[30] - x125*H[26] + 1;
   const float x130 = x74 + x75 + x76;
   const float x131 = x77 + x78 + x79;
   const float x132 = x80 + x81 + x82;
   const float x133 = -x130*H[33] - x131*H[29] - x132*H[25];
   const float x134 = -x130*H[34] - x131*H[30] - x132*H[26];
   const float x135 = -x130*H[36] - x131*H[32] - x132*H[28];
   const float x136 = -x130*H[35] - x131*H[31] - x132*H[27] + 1;
   const float x137 = x83 + x84 + x85;
   const float x138 = x86 + x87 + x88;
   const float x139 = x89 + x90 + x91;
   const float x140 = -x137*H[33] - x138*H[29] - x139*H[25];
   const float x141 = -x137*H[34] - x138*H[30] - x139*H[26];
   const float x142 = -x137*H[35] - x138*H[31] - x139*H[27];
   const float x143 = -x137*H[36] - x138*H[32] - x139*H[28] + 1;
   const float x144 = x93 + x95 + x97;
   const float x145 = x100 + x98 + x99;
   const float x146 = x101 + x102 + x103;
   const float x147 = -x144*H[33] - x145*H[29] - x146*H[25];
   const float x148 = -x144*H[34] - x145*H[30] - x146*H[26];
   const float x149 = -x144*H[35] - x145*H[31] - x146*H[27];
   const float x150 = -x144*H[36] - x145*H[32] - x146*H[28];
   const float x151 = x105 + x107 + x109;
   const float x152 = x110 + x111 + x112;
   const float x153 = x113 + x114 + x115;
   xnew[0] = x[0];
   xnew[1] = x[1];
   xnew[2] = x[2];
   xnew[3] = x[3];
   xnew[4] = x12*(-x42 - x44 - x46) + x51*(-x53 - x55 - x57) + x58*(-x60 - x62 - x64) + x[4];
   xnew[5] = x12*(-x65 - x66 - x67) + x51*(-x68 - x69 - x70) + x58*(-x71 - x72 - x73) + x[5];
   xnew[6] = x12*(-x74 - x75 - x76) + x51*(-x77 - x78 - x79) + x58*(-x80 - x81 - x82) + x[6];
   xnew[7] = x12*(-x83 - x84 - x85) + x51*(-x86 - x87 - x88) + x58*(-x89 - x90 - x91) + x[7];
   xnew[8] = x[8];
   xnew[9] = x[9];
   xnew[10] = x12*(-x93 - x95 - x97) + x51*(-x100 - x98 - x99) + x58*(-x101 - x102 - x103) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   xnew[18] = x[18];
   xnew[19] = x[19];
   xnew[20] = x12*(-x105 - x107 - x109) + x51*(-x110 - x111 - x112) + x58*(-x113 - x114 - x115) + x[20];
   xnew[21] = x[21];
   xnew[22] = x[22];
   Pnew[0] = P[0];
   Pnew[1] = P[1];
   Pnew[2] = P[2];
   Pnew[3] = P[3];
   Pnew[4] = P[4];
   Pnew[5] = P[5];
   Pnew[6] = P[6];
   Pnew[7] = P[7];
   Pnew[8] = P[8];
   Pnew[9] = P[9];
   Pnew[10] = P[10];
   Pnew[11] = P[11];
   Pnew[12] = P[12];
   Pnew[13] = P[13];
   Pnew[14] = P[14];
   Pnew[15] = P[15];
   Pnew[16] = P[16];
   Pnew[17] = P[17];
   Pnew[18] = P[18];
   Pnew[19] = P[19];
   Pnew[20] = P[20];
   Pnew[21] = P[21];
   Pnew[22] = P[22];
   Pnew[23] = P[23];
   Pnew[24] = P[24];
   Pnew[25] = P[25];
   Pnew[26] = P[26];
   Pnew[27] = P[27];
   Pnew[28] = P[28];
   Pnew[29] = P[29];
   Pnew[30] = P[30];
   Pnew[31] = P[31];
   Pnew[32] = P[32];
   Pnew[33] = P[33];
   Pnew[34] = P[34];
   Pnew[35] = P[35];
   Pnew[36] = P[36];
   Pnew[37] = P[37];
   Pnew[38] = P[38];
   Pnew[39] = P[39];
   Pnew[40] = P[40];
   Pnew[41] = P[41];
   Pnew[42] = P[42];
   Pnew[43] = P[43];
   Pnew[44] = P[44];
   Pnew[45] = P[45];
   Pnew[46] = P[46];
   Pnew[47] = P[47];
   Pnew[48] = P[48];
   Pnew[49] = P[49];
   Pnew[50] = x119*P[51] + x120*P[52] + x121*P[53] + x122*P[50];
   Pnew[51] = x119*P[63] + x120*P[64] + x121*P[65] + x122*P[51];
   Pnew[52] = x119*P[64] + x120*P[75] + x121*P[76] + x122*P[52];
   Pnew[53] = x119*P[65] + x120*P[76] + x121*P[86] + x122*P[53];
   Pnew[54] = P[54];
   Pnew[55] = P[55];
   Pnew[56] = x119*P[68] + x120*P[79] + x121*P[89] + x122*P[56];
   Pnew[57] = P[57];
   Pnew[58] = P[58];
   Pnew[59] = P[59];
   Pnew[60] = x119*P[72] + x120*P[83] + x121*P[93] + x122*P[60];
   Pnew[61] = P[61];
   Pnew[62] = P[62];
   Pnew[63] = x126*P[51] + x127*P[64] + x128*P[65] + x129*P[63];
   Pnew[64] = x126*P[52] + x127*P[75] + x128*P[76] + x129*P[64];
   Pnew[65] = x126*P[53] + x127*P[76] + x128*P[86] + x129*P[65];
   Pnew[66] = P[66];
   Pnew[67] = P[67];
   Pnew[68] = x126*P[56] + x127*P[79] + x128*P[89] + x129*P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = x126*P[60] + x127*P[83] + x128*P[93] + x129*P[72];
   Pnew[73] = P[73];
   Pnew[74] = P[74];
   Pnew[75] = x133*P[52] + x134*P[64] + x135*P[76] + x136*P[75];
   Pnew[76] = x133*P[53] + x134*P[65] + x135*P[86] + x136*P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = x133*P[56] + x134*P[68] + x135*P[89] + x136*P[79];
   Pnew[80] = P[80];
   Pnew[81] = P[81];
   Pnew[82] = P[82];
   Pnew[83] = x133*P[60] + x134*P[72] + x135*P[93] + x136*P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = x140*P[53] + x141*P[65] + x142*P[76] + x143*P[86];
   Pnew[87] = P[87];
   Pnew[88] = P[88];
   Pnew[89] = x140*P[56] + x141*P[68] + x142*P[79] + x143*P[89];
   Pnew[90] = P[90];
   Pnew[91] = P[91];
   Pnew[92] = P[92];
   Pnew[93] = x140*P[60] + x141*P[72] + x142*P[83] + x143*P[93];
   Pnew[94] = P[94];
   Pnew[95] = P[95];
   Pnew[96] = P[96];
   Pnew[97] = P[97];
   Pnew[98] = P[98];
   Pnew[99] = P[99];
   Pnew[100] = P[100];
   Pnew[101] = P[101];
   Pnew[102] = P[102];
   Pnew[103] = P[103];
   Pnew[104] = P[104];
   Pnew[105] = P[105];
   Pnew[106] = x147*P[56] + x148*P[68] + x149*P[79] + x150*P[89] + P[106];
   Pnew[107] = P[107];
   Pnew[108] = P[108];
   Pnew[109] = x147*P[60] + x148*P[72] + x149*P[83] + x150*P[93] + P[109];
   Pnew[110] = P[110];
   Pnew[111] = P[111];
   Pnew[112] = P[112];
   Pnew[113] = P[113];
   Pnew[114] = P[114];
   Pnew[115] = P[115];
   Pnew[116] = P[116];
   Pnew[117] = P[117];
   Pnew[118] = P[118];
   Pnew[119] = P[119];
   Pnew[120] = P[120];
   Pnew[121] = P[121];
   Pnew[122] = P[122];
   Pnew[123] = P[123];
   Pnew[124] = P[124];
   Pnew[125] = P[125];
   Pnew[126] = P[126];
   Pnew[127] = P[127];
   Pnew[128] = P[128];
   Pnew[129] = P[129];
   Pnew[130] = P[130];
   Pnew[131] = P[131];
   Pnew[132] = P[132];
   Pnew[133] = P[133];
   Pnew[134] = P[134];
   Pnew[135] = P[135];
   Pnew[136] = P[136];
   Pnew[137] = P[137];
   Pnew[138] = P[138];
   Pnew[139] = (-x151*H[33] - x152*H[29] - x153*H[25])*P[60] + (-x151*H[34] - x152*H[30] - x153*H[26])*P[72] + (-x151*H[35] - x152*H[31] - x153*H[27])*P[83] + (-x151*H[36] - x152*H[32] - x153*H[28])*P[93] + P[139];
   Pnew[140] = P[140];
   Pnew[141] = P[141];
   Pnew[142] = P[142];
   Pnew[143] = P[143];
}

void gyro_correction(const float * restrict x, const float * restrict P, const float * restrict gyro,
	const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = H[19]*P[96];
   const float x1 = H[20]*P[100];
   const float x2 = x0 + x1;
   const float x3 = x2*H[19];
   const float x4 = H[19]*P[100];
   const float x5 = H[20]*P[135];
   const float x6 = x4 + x5;
   const float x7 = x6*H[20];
   const float x8 = 1.0F/(x3 + x7 + R[4]);
   const float x9 = x8*(-gyro[0] + x[8] + x[18]);
   const float x10 = x9*H[19];
   const float x11 = H[21]*P[102];
   const float x12 = H[22]*P[104];
   const float x13 = x11 + x12;
   const float x14 = x13*H[21];
   const float x15 = H[21]*P[104];
   const float x16 = H[22]*P[137];
   const float x17 = x15 + x16;
   const float x18 = x17*H[22];
   const float x19 = 1.0F/(x14 + x18 + R[5]);
   const float x20 = x19*(-gyro[1] + x[9] + x[19]);
   const float x21 = x20*H[21];
   const float x22 = H[23]*P[106];
   const float x23 = H[24]*P[109];
   const float x24 = x22 + x23;
   const float x25 = x24*H[23];
   const float x26 = H[23]*P[109];
   const float x27 = H[24]*P[139];
   const float x28 = x26 + x27;
   const float x29 = x28*H[24];
   const float x30 = 1.0F/(x25 + x29 + R[6]);
   const float x31 = x30*(-gyro[2] + x[10] + x[20]);
   const float x32 = x31*H[23];
   const float x33 = H[19]*P[54];
   const float x34 = H[20]*P[58];
   const float x35 = x33 + x34;
   const float x36 = H[21]*P[55];
   const float x37 = H[22]*P[59];
   const float x38 = x36 + x37;
   const float x39 = H[23]*P[56];
   const float x40 = H[24]*P[60];
   const float x41 = x39 + x40;
   const float x42 = H[19]*P[66];
   const float x43 = H[20]*P[70];
   const float x44 = x42 + x43;
   const float x45 = H[21]*P[67];
   const float x46 = H[22]*P[71];
   const float x47 = x45 + x46;
   const float x48 = H[23]*P[68];
   const float x49 = H[24]*P[72];
   const float x50 = x48 + x49;
   const float x51 = H[19]*P[77];
   const float x52 = H[20]*P[81];
   const float x53 = x51 + x52;
   const float x54 = H[21]*P[78];
   const float x55 = H[22]*P[82];
   const float x56 = x54 + x55;
   const float x57 = H[23]*P[79];
   const float x58 = H[24]*P[83];
   const float x59 = x57 + x58;
   const float x60 = H[19]*P[87];
   const float x61 = H[20]*P[91];
   const float x62 = x60 + x61;
   const float x63 = H[21]*P[88];
   const float x64 = H[22]*P[92];
   const float x65 = x63 + x64;
   const float x66 = H[23]*P[89];
   const float x67 = H[24]*P[93];
   const float x68 = x66 + x67;
   const float x69 = H[19]*P[97];
   const float x70 = H[20]*P[113];
   const float x71 = x69 + x70;
   const float x72 = H[22]*P[117];
   const float x73 = H[23]*P[107];
   const float x74 = H[24]*P[121];
   const float x75 = x73 + x74;
   const float x76 = H[19]*P[99];
   const float x77 = H[20]*P[127];
   const float x78 = x76 + x77;
   const float x79 = H[21]*P[103];
   const float x80 = H[22]*P[130];
   const float x81 = x79 + x80;
   const float x82 = H[23]*P[108];
   const float x83 = H[24]*P[133];
   const float x84 = x82 + x83;
   const float x85 = H[19]*P[101];
   const float x86 = H[20]*P[136];
   const float x87 = x85 + x86;
   const float x88 = H[21]*P[105];
   const float x89 = H[22]*P[138];
   const float x90 = x88 + x89;
   const float x91 = H[23]*P[110];
   const float x92 = H[24]*P[140];
   const float x93 = x91 + x92;
   const float x94 = x8*(H[19]*H[19]);
   const float x95 = x19*(H[21]*H[21]);
   const float x96 = x30*(H[23]*H[23]);
   const float x97 = x94*P[8];
   const float x98 = x95*P[9];
   const float x99 = x96*P[10];
   const float x100 = x8*H[19];
   const float x101 = x100*P[8];
   const float x102 = x19*H[21];
   const float x103 = x102*P[9];
   const float x104 = x30*H[23];
   const float x105 = x104*P[10];
   const float x106 = x94*P[21];
   const float x107 = x95*P[22];
   const float x108 = x96*P[23];
   const float x109 = x100*P[21];
   const float x110 = x102*P[22];
   const float x111 = x104*P[23];
   const float x112 = x94*P[33];
   const float x113 = x95*P[34];
   const float x114 = x96*P[35];
   const float x115 = x100*P[33];
   const float x116 = x102*P[34];
   const float x117 = x104*P[35];
   const float x118 = x94*P[44];
   const float x119 = x95*P[45];
   const float x120 = x96*P[46];
   const float x121 = x100*P[44];
   const float x122 = x102*P[45];
   const float x123 = x104*P[46];
   const float x124 = x35*x8;
   const float x125 = x19*x38;
   const float x126 = x30*x41;
   const float x127 = x44*x8;
   const float x128 = x19*x47;
   const float x129 = x30*x50;
   const float x130 = x53*x8;
   const float x131 = x19*x56;
   const float x132 = x30*x59;
   const float x133 = x62*x8;
   const float x134 = x19*x65;
   const float x135 = x30*x68;
   const float x136 = x2*x8;
   const float x137 = -x3*x8 + 1;
   const float x138 = x13*x19;
   const float x139 = -x14*x19 + 1;
   const float x140 = x24*x30;
   const float x141 = -x25*x30 + 1;
   const float x142 = x71*x8;
   const float x143 = x19*(H[22]*H[22]);
   const float x144 = x143*P[117];
   const float x145 = x19*x72;
   const float x146 = x30*x75;
   const float x147 = x78*x8;
   const float x148 = x19*x81;
   const float x149 = x30*x84;
   const float x150 = x6*x8;
   const float x151 = -x7*x8 + 1;
   const float x152 = x17*x19;
   const float x153 = -x18*x19 + 1;
   const float x154 = x28*x30;
   const float x155 = -x29*x30 + 1;
   const float x156 = x8*x87;
   const float x157 = x19*x90;
   const float x158 = x30*x93;
   xnew[0] = -x10*P[8] - x21*P[9] - x32*P[10] + x[0];
   xnew[1] = -x10*P[21] - x21*P[22] - x32*P[23] + x[1];
   xnew[2] = -x10*P[33] - x21*P[34] - x32*P[35] + x[2];
   xnew[3] = -x10*P[44] - x21*P[45] - x32*P[46] + x[3];
   xnew[4] = -x20*x38 - x31*x41 - x35*x9 + x[4];
   xnew[5] = -x20*x47 - x31*x50 - x44*x9 + x[5];
   xnew[6] = -x20*x56 - x31*x59 - x53*x9 + x[6];
   xnew[7] = -x20*x65 - x31*x68 - x62*x9 + x[7];
   xnew[8] = -x2*x9 + x[8];
   xnew[9] = -x13*x20 + x[9];
   xnew[10] = -x24*x31 + x[10];
   xnew[11] = -x71*x9 + x[11];
   xnew[12] = -x10*P[98] - x20*x72 + x[12];
   xnew[13] = -x31*x75 + x[13];
   xnew[14] = x[14];
   xnew[15] = -x78*x9 + x[15];
   xnew[16] = -x20*x81 + x[16];
   xnew[17] = -x31*x84 + x[17];
   xnew[18] = -x6*x9 + x[18];
   xnew[19] = -x17*x20 + x[19];
   xnew[20] = -x28*x31 + x[20];
   xnew[21] = x[21];
   xnew[22] = -x20*x90 - x31*x93 - x87*x9 + x[22];
   Pnew[0] = -x94*(P[8]*P[8]) - x95*(P[9]*P[9]) - x96*(P[10]*P[10]) + P[0];
   Pnew[1] = -x97*P[21] - x98*P[22] - x99*P[23] + P[1];
   Pnew[2] = -x97*P[33] - x98*P[34] - x99*P[35] + P[2];
   Pnew[3] = -x97*P[44] - x98*P[45] - x99*P[46] + P[3];
   Pnew[4] = -x101*x34 - x103*x37 - x105*x40 - x97*P[54] - x98*P[55] - x99*P[56] + P[4];
   Pnew[5] = -x101*x43 - x103*x46 - x105*x49 - x97*P[66] - x98*P[67] - x99*P[68] + P[5];
   Pnew[6] = -x101*x52 - x103*x55 - x105*x58 - x97*P[77] - x98*P[78] - x99*P[79] + P[6];
   Pnew[7] = -x101*x61 - x103*x64 - x105*x67 - x97*P[87] - x98*P[88] - x99*P[89] + P[7];
   Pnew[8] = -x1*x101 - x97*P[96] + P[8];
   Pnew[9] = -x103*x12 - x98*P[102] + P[9];
   Pnew[10] = -x105*x23 - x99*P[106] + P[10];
   Pnew[11] = P[11];
   Pnew[12] = P[12];
   Pnew[13] = -x101*x86 - x103*x89 - x105*x92 - x97*P[101] - x98*P[105] - x99*P[110] + P[13];
   Pnew[14] = -x94*(P[21]*P[21]) - x95*(P[22]*P[22]) - x96*(P[23]*P[23]) + P[14];
   Pnew[15] = -x106*P[33] - x107*P[34] - x108*P[35] + P[15];
   Pnew[16] = -x106*P[44] - x107*P[45] - x108*P[46] + P[16];
   Pnew[17] = -x106*P[54] - x107*P[55] - x108*P[56] - x109*x34 - x110*x37 - x111*x40 + P[17];
   Pnew[18] = -x106*P[66] - x107*P[67] - x108*P[68] - x109*x43 - x110*x46 - x111*x49 + P[18];
   Pnew[19] = -x106*P[77] - x107*P[78] - x108*P[79] - x109*x52 - x110*x55 - x111*x58 + P[19];
   Pnew[20] = -x106*P[87] - x107*P[88] - x108*P[89] - x109*x61 - x110*x64 - x111*x67 + P[20];
   Pnew[21] = -x1*x109 - x106*P[96] + P[21];
   Pnew[22] = -x107*P[102] - x110*x12 + P[22];
   Pnew[23] = -x108*P[106] - x111*x23 + P[23];
   Pnew[24] = P[24];
   Pnew[25] = P[25];
   Pnew[26] = -x106*P[101] - x107*P[105] - x108*P[110] - x109*x86 - x110*x89 - x111*x92 + P[26];
   Pnew[27] = -x94*(P[33]*P[33]) - x95*(P[34]*P[34]) - x96*(P[35]*P[35]) + P[27];
   Pnew[28] = -x112*P[44] - x113*P[45] - x114*P[46] + P[28];
   Pnew[29] = -x112*P[54] - x113*P[55] - x114*P[56] - x115*x34 - x116*x37 - x117*x40 + P[29];
   Pnew[30] = -x112*P[66] - x113*P[67] - x114*P[68] - x115*x43 - x116*x46 - x117*x49 + P[30];
   Pnew[31] = -x112*P[77] - x113*P[78] - x114*P[79] - x115*x52 - x116*x55 - x117*x58 + P[31];
   Pnew[32] = -x112*P[87] - x113*P[88] - x114*P[89] - x115*x61 - x116*x64 - x117*x67 + P[32];
   Pnew[33] = -x1*x115 - x112*P[96] + P[33];
   Pnew[34] = -x113*P[102] - x116*x12 + P[34];
   Pnew[35] = -x114*P[106] - x117*x23 + P[35];
   Pnew[36] = P[36];
   Pnew[37] = P[37];
   Pnew[38] = -x112*P[101] - x113*P[105] - x114*P[110] - x115*x86 - x116*x89 - x117*x92 + P[38];
   Pnew[39] = -x94*(P[44]*P[44]) - x95*(P[45]*P[45]) - x96*(P[46]*P[46]) + P[39];
   Pnew[40] = -x118*P[54] - x119*P[55] - x120*P[56] - x121*x34 - x122*x37 - x123*x40 + P[40];
   Pnew[41] = -x118*P[66] - x119*P[67] - x120*P[68] - x121*x43 - x122*x46 - x123*x49 + P[41];
   Pnew[42] = -x118*P[77] - x119*P[78] - x120*P[79] - x121*x52 - x122*x55 - x123*x58 + P[42];
   Pnew[43] = -x118*P[87] - x119*P[88] - x120*P[89] - x121*x61 - x122*x64 - x123*x67 + P[43];
   Pnew[44] = -x1*x121 - x118*P[96] + P[44];
   Pnew[45] = -x119*P[102] - x12*x122 + P[45];
   Pnew[46] = -x120*P[106] - x123*x23 + P[46];
   Pnew[47] = P[47];
   Pnew[48] = P[48];
   Pnew[49] = -x118*P[101] - x119*P[105] - x120*P[110] - x121*x86 - x122*x89 - x123*x92 + P[49];
   Pnew[50] = -x124*x33 - x124*x34 - x125*x36 - x125*x37 - x126*x39 - x126*x40 + P[50];
   Pnew[51] = -x124*x42 - x124*x43 - x125*x45 - x125*x46 - x126*x48 - x126*x49 + P[51];
   Pnew[52] = -x124*x51 - x124*x52 - x125*x54 - x125*x55 - x126*x57 - x126*x58 + P[52];
   Pnew[53] = -x124*x60 - x124*x61 - x125*x63 - x125*x64 - x126*x66 - x126*x67 + P[53];
   Pnew[54] = -x0*x124 - x1*x124 + P[54];
   Pnew[55] = -x11*x125 - x12*x125 + P[55];
   Pnew[56] = -x126*x22 - x126*x23 + P[56];
   Pnew[57] = P[57];
   Pnew[58] = -x124*x4 - x124*x5 + P[58];
   Pnew[59] = -x125*x15 - x125*x16 + P[59];
   Pnew[60] = -x126*x26 - x126*x27 + P[60];
   Pnew[61] = P[61];
   Pnew[62] = -x124*x85 - x124*x86 - x125*x88 - x125*x89 - x126*x91 - x126*x92 + P[62];
   Pnew[63] = -x127*x42 - x127*x43 - x128*x45 - x128*x46 - x129*x48 - x129*x49 + P[63];
   Pnew[64] = -x127*x51 - x127*x52 - x128*x54 - x128*x55 - x129*x57 - x129*x58 + P[64];
   Pnew[65] = -x127*x60 - x127*x61 - x128*x63 - x128*x64 - x129*x66 - x129*x67 + P[65];
   Pnew[66] = -x0*x127 - x1*x127 + P[66];
   Pnew[67] = -x11*x128 - x12*x128 + P[67];
   Pnew[68] = -x129*x22 - x129*x23 + P[68];
   Pnew[69] = P[69];
   Pnew[70] = -x127*x4 - x127*x5 + P[70];
   Pnew[71] = -x128*x15 - x128*x16 + P[71];
   Pnew[72] = -x129*x26 - x129*x27 + P[72];
   Pnew[73] = P[73];
   Pnew[74] = -x127*x85 - x127*x86 - x128*x88 - x128*x89 - x129*x91 - x129*x92 + P[74];
   Pnew[75] = -x130*x51 - x130*x52 - x131*x54 - x131*x55 - x132*x57 - x132*x58 + P[75];
   Pnew[76] = -x130*x60 - x130*x61 - x131*x63 - x131*x64 - x132*x66 - x132*x67 + P[76];
   Pnew[77] = -x0*x130 - x1*x130 + P[77];
   Pnew[78] = -x11*x131 - x12*x131 + P[78];
   Pnew[79] = -x132*x22 - x132*x23 + P[79];
   Pnew[80] = P[80];
   Pnew[81] = -x130*x4 - x130*x5 + P[81];
   Pnew[82] = -x131*x15 - x131*x16 + P[82];
   Pnew[83] = -x132*x26 - x132*x27 + P[83];
   Pnew[84] = P[84];
   Pnew[85] = -x130*x85 - x130*x86 - x131*x88 - x131*x89 - x132*x91 - x132*x92 + P[85];
   Pnew[86] = -x133*x60 - x133*x61 - x134*x63 - x134*x64 - x135*x66 - x135*x67 + P[86];
   Pnew[87] = -x0*x133 - x1*x133 + P[87];
   Pnew[88] = -x11*x134 - x12*x134 + P[88];
   Pnew[89] = -x135*x22 - x135*x23 + P[89];
   Pnew[90] = P[90];
   Pnew[91] = -x133*x4 - x133*x5 + P[91];
   Pnew[92] = -x134*x15 - x134*x16 + P[92];
   Pnew[93] = -x135*x26 - x135*x27 + P[93];
   Pnew[94] = P[94];
   Pnew[95] = -x133*x85 - x133*x86 - x134*x88 - x134*x89 - x135*x91 - x135*x92 + P[95];
   Pnew[96] = -x1*x136 + x137*P[96];
   Pnew[97] = -x136*x70 + x137*P[97];
   Pnew[98] = x137*P[98];
   Pnew[99] = -x136*x77 + x137*P[99];
   Pnew[100] = -x136*x5 + x137*P[100];
   Pnew[101] = -x136*x86 + x137*P[101];
   Pnew[102] = -x12*x138 + x139*P[102];
   Pnew[103] = -x138*x80 + x139*P[103];
   Pnew[104] = -x138*x16 + x139*P[104];
   Pnew[105] = -x138*x89 + x139*P[105];
   Pnew[106] = -x140*x23 + x141*P[106];
   Pnew[107] = -x140*x74 + x141*P[107];
   Pnew[108] = -x140*x83 + x141*P[108];
   Pnew[109] = -x140*x27 + x141*P[109];
   Pnew[110] = -x140*x92 + x141*P[110];
   Pnew[111] = -x142*x69 - x142*x70 + P[111];
   Pnew[112] = -x142*x76 - x142*x77 + P[112];
   Pnew[113] = -x142*x4 - x142*x5 + P[113];
   Pnew[114] = -x142*x85 - x142*x86 + P[114];
   Pnew[115] = -x143*(P[117]*P[117]) - x94*(P[98]*P[98]) + P[115];
   Pnew[116] = -x144*P[130] - x145*x79 + P[116];
   Pnew[117] = -x102*x12*P[117] - x144*P[137] + P[117];
   Pnew[118] = -x100*x86*P[98] - x144*P[138] - x145*x88 - x94*P[98]*P[101] + P[118];
   Pnew[119] = -x146*x73 - x146*x74 + P[119];
   Pnew[120] = -x146*x82 - x146*x83 + P[120];
   Pnew[121] = -x146*x26 - x146*x27 + P[121];
   Pnew[122] = -x146*x91 - x146*x92 + P[122];
   Pnew[123] = P[123];
   Pnew[124] = P[124];
   Pnew[125] = P[125];
   Pnew[126] = -x147*x76 - x147*x77 + P[126];
   Pnew[127] = -x147*x4 - x147*x5 + P[127];
   Pnew[128] = -x147*x85 - x147*x86 + P[128];
   Pnew[129] = -x148*x79 - x148*x80 + P[129];
   Pnew[130] = -x148*x15 - x148*x16 + P[130];
   Pnew[131] = -x148*x88 - x148*x89 + P[131];
   Pnew[132] = -x149*x82 - x149*x83 + P[132];
   Pnew[133] = -x149*x26 - x149*x27 + P[133];
   Pnew[134] = -x149*x91 - x149*x92 + P[134];
   Pnew[135] = -x150*x4 + x151*P[135];
   Pnew[136] = -x150*x85 + x151*P[136];
   Pnew[137] = -x15*x152 + x153*P[137];
   Pnew[138] = -x152*x88 + x153*P[138];
   Pnew[139] = -x154*x26 + x155*P[139];
   Pnew[140] = -x154*x91 + x155*P[140];
   Pnew[141] = P[141];
   Pnew[142] = P[142];
   Pnew[143] = -x156*x85 - x156*x86 - x157*x88 - x157*x89 - x158*x91 - x158*x92 + P[143];
}

void accel_correction(const float *restrict x, const float *restrict P, const float *restrict accel,
        const float *restrict H, const float *restrict R, const float *restrict param, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = expf(x[22]);
   const float x1 = 2*x[3];
   const float x2 = x[4]*x[7];
   const float x3 = x[5]*x[6];
   const float x4 = (x[6]*x[6]);
   const float x5 = (x[5]*x[5]);
   const float x6 = (x[4]*x[4]) - (x[7]*x[7]);
   const float x7 = x0*(-x1*(x[4]*x[5] + x[6]*x[7]) + 2*(x2 - x3)*x[1] - (x4 - x5 + x6)*x[2]) - accel[1];
   const float x8 = H[1]*P[24] + H[2]*P[36] + H[3]*P[47] + H[4]*P[57] + H[5]*P[69] + H[6]*P[80] + H[7]*P[90] + H[8]*P[125];
   const float x9 = H[1]*P[25] + H[2]*P[37] + H[3]*P[48] + H[4]*P[61] + H[5]*P[73] + H[6]*P[84] + H[7]*P[94] + H[8]*P[142];
   const float x10 = x8*H[17] + x9*H[18];
   const float x11 = H[17]*P[24];
   const float x12 = H[18]*P[25];
   const float x13 = x11 + x12;
   const float x14 = H[17]*P[36];
   const float x15 = H[18]*P[37];
   const float x16 = x14 + x15;
   const float x17 = H[17]*P[47];
   const float x18 = H[18]*P[48];
   const float x19 = x17 + x18;
   const float x20 = H[17]*P[57];
   const float x21 = H[18]*P[61];
   const float x22 = x20 + x21;
   const float x23 = H[17]*P[69];
   const float x24 = H[18]*P[73];
   const float x25 = x23 + x24;
   const float x26 = H[17]*P[80];
   const float x27 = H[18]*P[84];
   const float x28 = x26 + x27;
   const float x29 = H[17]*P[90];
   const float x30 = H[18]*P[94];
   const float x31 = x29 + x30;
   const float x32 = H[17]*P[125];
   const float x33 = H[18]*P[142];
   const float x34 = x32 + x33;
   const float x35 = x13*H[9] + x16*H[10] + x19*H[11] + x22*H[12] + x25*H[13] + x28*H[14] + x31*H[15] + x34*H[16];
   const float x36 = H[17]*P[123];
   const float x37 = H[18]*P[124];
   const float x38 = x36 + x37;
   const float x39 = H[17]*P[124];
   const float x40 = H[18]*P[141];
   const float x41 = x39 + x40;
   const float x42 = x38*H[17] + x41*H[18] + R[3];
   const float x43 = H[1]*P[14] + H[2]*P[15] + H[3]*P[16] + H[4]*P[17] + H[5]*P[18] + H[6]*P[19] + H[7]*P[20] + H[8]*P[26];
   const float x44 = H[1]*P[15] + H[2]*P[27] + H[3]*P[28] + H[4]*P[29] + H[5]*P[30] + H[6]*P[31] + H[7]*P[32] + H[8]*P[38];
   const float x45 = H[1]*P[16] + H[2]*P[28] + H[3]*P[39] + H[4]*P[40] + H[5]*P[41] + H[6]*P[42] + H[7]*P[43] + H[8]*P[49];
   const float x46 = H[1]*P[17] + H[2]*P[29] + H[3]*P[40] + H[4]*P[50] + H[5]*P[51] + H[6]*P[52] + H[7]*P[53] + H[8]*P[62];
   const float x47 = H[1]*P[18] + H[2]*P[30] + H[3]*P[41] + H[4]*P[51] + H[5]*P[63] + H[6]*P[64] + H[7]*P[65] + H[8]*P[74];
   const float x48 = H[1]*P[19] + H[2]*P[31] + H[3]*P[42] + H[4]*P[52] + H[5]*P[64] + H[6]*P[75] + H[7]*P[76] + H[8]*P[85];
   const float x49 = H[1]*P[20] + H[2]*P[32] + H[3]*P[43] + H[4]*P[53] + H[5]*P[65] + H[6]*P[76] + H[7]*P[86] + H[8]*P[95];
   const float x50 = H[1]*P[26] + H[2]*P[38] + H[3]*P[49] + H[4]*P[62] + H[5]*P[74] + H[6]*P[85] + H[7]*P[95] + H[8]*P[143];
   const float x51 = x43*H[9] + x44*H[10] + x45*H[11] + x46*H[12] + x47*H[13] + x48*H[14] + x49*H[15] + x50*H[16];
   const float x52 = x10*x35 - x42*x51;
   const float x53 = H[9]*P[14] + H[10]*P[15] + H[11]*P[16] + H[12]*P[17] + H[13]*P[18] + H[14]*P[19] + H[15]*P[20] + H[16]*P[26];
   const float x54 = H[9]*P[15] + H[10]*P[27] + H[11]*P[28] + H[12]*P[29] + H[13]*P[30] + H[14]*P[31] + H[15]*P[32] + H[16]*P[38];
   const float x55 = H[9]*P[16] + H[10]*P[28] + H[11]*P[39] + H[12]*P[40] + H[13]*P[41] + H[14]*P[42] + H[15]*P[43] + H[16]*P[49];
   const float x56 = H[9]*P[17] + H[10]*P[29] + H[11]*P[40] + H[12]*P[50] + H[13]*P[51] + H[14]*P[52] + H[15]*P[53] + H[16]*P[62];
   const float x57 = H[9]*P[18] + H[10]*P[30] + H[11]*P[41] + H[12]*P[51] + H[13]*P[63] + H[14]*P[64] + H[15]*P[65] + H[16]*P[74];
   const float x58 = H[9]*P[19] + H[10]*P[31] + H[11]*P[42] + H[12]*P[52] + H[13]*P[64] + H[14]*P[75] + H[15]*P[76] + H[16]*P[85];
   const float x59 = H[9]*P[20] + H[10]*P[32] + H[11]*P[43] + H[12]*P[53] + H[13]*P[65] + H[14]*P[76] + H[15]*P[86] + H[16]*P[95];
   const float x60 = H[9]*P[26] + H[10]*P[38] + H[11]*P[49] + H[12]*P[62] + H[13]*P[74] + H[14]*P[85] + H[15]*P[95] + H[16]*P[143];
   const float x61 = x53*H[1] + x54*H[2] + x55*H[3] + x56*H[4] + x57*H[5] + x58*H[6] + x59*H[7] + x60*H[8];
   const float x62 = x35*x61;
   const float x63 = H[9]*P[24] + H[10]*P[36] + H[11]*P[47] + H[12]*P[57] + H[13]*P[69] + H[14]*P[80] + H[15]*P[90] + H[16]*P[125];
   const float x64 = H[9]*P[25] + H[10]*P[37] + H[11]*P[48] + H[12]*P[61] + H[13]*P[73] + H[14]*P[84] + H[15]*P[94] + H[16]*P[142];
   const float x65 = x63*H[17] + x64*H[18];
   const float x66 = x13*H[1] + x16*H[2] + x19*H[3] + x22*H[4] + x25*H[5] + x28*H[6] + x31*H[7] + x34*H[8];
   const float x67 = x65*x66;
   const float x68 = x53*H[9] + x54*H[10] + x55*H[11] + x56*H[12] + x57*H[13] + x58*H[14] + x59*H[15] + x60*H[16] + R[2];
   const float x69 = x66*x68;
   const float x70 = x43*H[1] + x44*H[2] + x45*H[3] + x46*H[4] + x47*H[5] + x48*H[6] + x49*H[7] + x50*H[8] + R[1];
   const float x71 = x35*x65;
   const float x72 = x42*x61;
   const float x73 = x42*x68;
   const float x74 = 1.0F/(x10*x62 - x10*x69 + x51*x67 - x51*x72 - x70*x71 + x70*x73);
   const float x75 = x74*(H[1]*P[1] + H[2]*P[2] + H[3]*P[3] + H[4]*P[4] + H[5]*P[5] + H[6]*P[6] + H[7]*P[7] + H[8]*P[13]);
   const float x76 = x52*x75;
   const float x77 = -x10*x66 + x42*x70;
   const float x78 = x74*(H[9]*P[1] + H[10]*P[2] + H[11]*P[3] + H[12]*P[4] + H[13]*P[5] + H[14]*P[6] + H[15]*P[7] + H[16]*P[13]);
   const float x79 = x77*x78;
   const float x80 = -x35*x70 + x51*x66;
   const float x81 = H[17]*P[11];
   const float x82 = H[18]*P[12];
   const float x83 = x74*(x81 + x82);
   const float x84 = x80*x83;
   const float x85 = x0*(x1*(x[4]*x[6] - x[5]*x[7]) - 2*(x2 + x3)*x[2] - (-x4 + x5 + x6)*x[1]) - accel[0];
   const float x86 = x67 - x72;
   const float x87 = x78*x86;
   const float x88 = -x71 + x73;
   const float x89 = x75*x88;
   const float x90 = x62 - x69;
   const float x91 = x83*x90;
   const float x92 = -expf(x[21])*param[0]*x[14] - accel[2];
   const float x93 = x10*x61 - x65*x70;
   const float x94 = x78*x93;
   const float x95 = -x10*x68 + x51*x65;
   const float x96 = x75*x95;
   const float x97 = -x51*x61 + x68*x70;
   const float x98 = x83*x97;
   const float x99 = x43*x74;
   const float x100 = x52*x99;
   const float x101 = x53*x74;
   const float x102 = x101*x77;
   const float x103 = x13*x74;
   const float x104 = x103*x80;
   const float x105 = x101*x86;
   const float x106 = x88*x99;
   const float x107 = x103*x90;
   const float x108 = x101*x93;
   const float x109 = x95*x99;
   const float x110 = x103*x97;
   const float x111 = x44*x74;
   const float x112 = x111*x52;
   const float x113 = x54*x74;
   const float x114 = x113*x77;
   const float x115 = x16*x74;
   const float x116 = x115*x80;
   const float x117 = x113*x86;
   const float x118 = x111*x88;
   const float x119 = x115*x90;
   const float x120 = x113*x93;
   const float x121 = x111*x95;
   const float x122 = x115*x97;
   const float x123 = x45*x74;
   const float x124 = x123*x52;
   const float x125 = x55*x74;
   const float x126 = x125*x77;
   const float x127 = x19*x74;
   const float x128 = x127*x80;
   const float x129 = x125*x86;
   const float x130 = x123*x88;
   const float x131 = x127*x90;
   const float x132 = x125*x93;
   const float x133 = x123*x95;
   const float x134 = x127*x97;
   const float x135 = x46*x74;
   const float x136 = x135*x52;
   const float x137 = x56*x74;
   const float x138 = x137*x77;
   const float x139 = x22*x74;
   const float x140 = x139*x80;
   const float x141 = x137*x86;
   const float x142 = x135*x88;
   const float x143 = x139*x90;
   const float x144 = x137*x93;
   const float x145 = x135*x95;
   const float x146 = x139*x97;
   const float x147 = x47*x74;
   const float x148 = x147*x52;
   const float x149 = x57*x74;
   const float x150 = x149*x77;
   const float x151 = x25*x74;
   const float x152 = x151*x80;
   const float x153 = x149*x86;
   const float x154 = x147*x88;
   const float x155 = x151*x90;
   const float x156 = x149*x93;
   const float x157 = x147*x95;
   const float x158 = x151*x97;
   const float x159 = x48*x74;
   const float x160 = x159*x52;
   const float x161 = x58*x74;
   const float x162 = x161*x77;
   const float x163 = x28*x74;
   const float x164 = x163*x80;
   const float x165 = x161*x86;
   const float x166 = x159*x88;
   const float x167 = x163*x90;
   const float x168 = x161*x93;
   const float x169 = x159*x95;
   const float x170 = x163*x97;
   const float x171 = x49*x74;
   const float x172 = x171*x52;
   const float x173 = x59*x74;
   const float x174 = x173*x77;
   const float x175 = x31*x74;
   const float x176 = x175*x80;
   const float x177 = x173*x86;
   const float x178 = x171*x88;
   const float x179 = x175*x90;
   const float x180 = x173*x93;
   const float x181 = x171*x95;
   const float x182 = x175*x97;
   const float x183 = x74*(H[1]*P[21] + H[2]*P[33] + H[3]*P[44] + H[4]*P[54] + H[5]*P[66] + H[6]*P[77] + H[7]*P[87] + H[8]*P[101]);
   const float x184 = x183*x52;
   const float x185 = x74*(H[9]*P[21] + H[10]*P[33] + H[11]*P[44] + H[12]*P[54] + H[13]*P[66] + H[14]*P[77] + H[15]*P[87] + H[16]*P[101]);
   const float x186 = x185*x77;
   const float x187 = x185*x86;
   const float x188 = x183*x88;
   const float x189 = x185*x93;
   const float x190 = x183*x95;
   const float x191 = x74*(H[1]*P[22] + H[2]*P[34] + H[3]*P[45] + H[4]*P[55] + H[5]*P[67] + H[6]*P[78] + H[7]*P[88] + H[8]*P[105]);
   const float x192 = x191*x52;
   const float x193 = x74*(H[9]*P[22] + H[10]*P[34] + H[11]*P[45] + H[12]*P[55] + H[13]*P[67] + H[14]*P[78] + H[15]*P[88] + H[16]*P[105]);
   const float x194 = x193*x77;
   const float x195 = x193*x86;
   const float x196 = x191*x88;
   const float x197 = x193*x93;
   const float x198 = x191*x95;
   const float x199 = x74*(H[1]*P[23] + H[2]*P[35] + H[3]*P[46] + H[4]*P[56] + H[5]*P[68] + H[6]*P[79] + H[7]*P[89] + H[8]*P[110]);
   const float x200 = x199*x52;
   const float x201 = x74*(H[9]*P[23] + H[10]*P[35] + H[11]*P[46] + H[12]*P[56] + H[13]*P[68] + H[14]*P[79] + H[15]*P[89] + H[16]*P[110]);
   const float x202 = x201*x77;
   const float x203 = x201*x86;
   const float x204 = x199*x88;
   const float x205 = x201*x93;
   const float x206 = x199*x95;
   const float x207 = x74*P[114];
   const float x208 = x207*H[8];
   const float x209 = x208*x52;
   const float x210 = x207*H[16];
   const float x211 = x210*x77;
   const float x212 = x210*x86;
   const float x213 = x208*x88;
   const float x214 = x210*x93;
   const float x215 = x208*x95;
   const float x216 = x74*P[118];
   const float x217 = x216*H[8];
   const float x218 = x217*x52;
   const float x219 = x216*H[16];
   const float x220 = x219*x77;
   const float x221 = x219*x86;
   const float x222 = x217*x88;
   const float x223 = x219*x93;
   const float x224 = x217*x95;
   const float x225 = x74*P[122];
   const float x226 = x225*H[8];
   const float x227 = x226*x52;
   const float x228 = x225*H[16];
   const float x229 = x228*x77;
   const float x230 = x228*x86;
   const float x231 = x226*x88;
   const float x232 = x228*x93;
   const float x233 = x226*x95;
   const float x234 = x74*x8;
   const float x235 = x234*x52;
   const float x236 = x63*x74;
   const float x237 = x236*x77;
   const float x238 = x38*x74;
   const float x239 = x238*x80;
   const float x240 = x236*x86;
   const float x241 = x234*x88;
   const float x242 = x238*x90;
   const float x243 = x236*x93;
   const float x244 = x234*x95;
   const float x245 = x238*x97;
   const float x246 = x74*P[128];
   const float x247 = x246*H[8];
   const float x248 = x247*x52;
   const float x249 = x246*H[16];
   const float x250 = x249*x77;
   const float x251 = x249*x86;
   const float x252 = x247*x88;
   const float x253 = x249*x93;
   const float x254 = x247*x95;
   const float x255 = x74*P[131];
   const float x256 = x255*H[8];
   const float x257 = x256*x52;
   const float x258 = x255*H[16];
   const float x259 = x258*x77;
   const float x260 = x258*x86;
   const float x261 = x256*x88;
   const float x262 = x258*x93;
   const float x263 = x256*x95;
   const float x264 = x74*P[134];
   const float x265 = x264*H[8];
   const float x266 = x265*x52;
   const float x267 = x264*H[16];
   const float x268 = x267*x77;
   const float x269 = x267*x86;
   const float x270 = x265*x88;
   const float x271 = x267*x93;
   const float x272 = x265*x95;
   const float x273 = x74*(H[4]*P[58] + H[5]*P[70] + H[6]*P[81] + H[7]*P[91] + H[8]*P[136]);
   const float x274 = x273*x52;
   const float x275 = x74*(H[12]*P[58] + H[13]*P[70] + H[14]*P[81] + H[15]*P[91] + H[16]*P[136]);
   const float x276 = x275*x77;
   const float x277 = x275*x86;
   const float x278 = x273*x88;
   const float x279 = x275*x93;
   const float x280 = x273*x95;
   const float x281 = x74*(H[4]*P[59] + H[5]*P[71] + H[6]*P[82] + H[7]*P[92] + H[8]*P[138]);
   const float x282 = x281*x52;
   const float x283 = x74*(H[12]*P[59] + H[13]*P[71] + H[14]*P[82] + H[15]*P[92] + H[16]*P[138]);
   const float x284 = x283*x77;
   const float x285 = x283*x86;
   const float x286 = x281*x88;
   const float x287 = x283*x93;
   const float x288 = x281*x95;
   const float x289 = x74*(H[4]*P[60] + H[5]*P[72] + H[6]*P[83] + H[7]*P[93] + H[8]*P[140]);
   const float x290 = x289*x52;
   const float x291 = x74*(H[12]*P[60] + H[13]*P[72] + H[14]*P[83] + H[15]*P[93] + H[16]*P[140]);
   const float x292 = x291*x77;
   const float x293 = x291*x86;
   const float x294 = x289*x88;
   const float x295 = x291*x93;
   const float x296 = x289*x95;
   const float x297 = x74*x9;
   const float x298 = x297*x52;
   const float x299 = x64*x74;
   const float x300 = x299*x77;
   const float x301 = x41*x74;
   const float x302 = x301*x80;
   const float x303 = x299*x86;
   const float x304 = x297*x88;
   const float x305 = x301*x90;
   const float x306 = x299*x93;
   const float x307 = x297*x95;
   const float x308 = x301*x97;
   const float x309 = x50*x74;
   const float x310 = x309*x52;
   const float x311 = x60*x74;
   const float x312 = x311*x77;
   const float x313 = x34*x74;
   const float x314 = x313*x80;
   const float x315 = x311*x86;
   const float x316 = x309*x88;
   const float x317 = x313*x90;
   const float x318 = x311*x93;
   const float x319 = x309*x95;
   const float x320 = x313*x97;
   const float x321 = x94 + x96 + x98;
   const float x322 = x76 + x79 + x84;
   const float x323 = x87 + x89 + x91;
   const float x324 = -x322*H[9] - x323*H[1];
   const float x325 = -x322*H[10] - x323*H[2];
   const float x326 = -x322*H[11] - x323*H[3];
   const float x327 = -x322*H[12] - x323*H[4];
   const float x328 = -x322*H[13] - x323*H[5];
   const float x329 = -x322*H[14] - x323*H[6];
   const float x330 = -x322*H[15] - x323*H[7];
   const float x331 = -x322*H[16] - x323*H[8];
   const float x332 = x108 + x109 + x110;
   const float x333 = x100 + x102 + x104;
   const float x334 = x105 + x106 + x107;
   const float x335 = -x333*H[10] - x334*H[2];
   const float x336 = -x333*H[11] - x334*H[3];
   const float x337 = -x333*H[12] - x334*H[4];
   const float x338 = -x333*H[13] - x334*H[5];
   const float x339 = -x333*H[14] - x334*H[6];
   const float x340 = -x333*H[15] - x334*H[7];
   const float x341 = -x333*H[16] - x334*H[8];
   const float x342 = -x333*H[9] - x334*H[1] + 1;
   const float x343 = x120 + x121 + x122;
   const float x344 = x112 + x114 + x116;
   const float x345 = x117 + x118 + x119;
   const float x346 = -x344*H[9] - x345*H[1];
   const float x347 = -x344*H[11] - x345*H[3];
   const float x348 = -x344*H[12] - x345*H[4];
   const float x349 = -x344*H[13] - x345*H[5];
   const float x350 = -x344*H[14] - x345*H[6];
   const float x351 = -x344*H[15] - x345*H[7];
   const float x352 = -x344*H[16] - x345*H[8];
   const float x353 = -x344*H[10] - x345*H[2] + 1;
   const float x354 = x132 + x133 + x134;
   const float x355 = x124 + x126 + x128;
   const float x356 = x129 + x130 + x131;
   const float x357 = -x355*H[9] - x356*H[1];
   const float x358 = -x355*H[10] - x356*H[2];
   const float x359 = -x355*H[12] - x356*H[4];
   const float x360 = -x355*H[13] - x356*H[5];
   const float x361 = -x355*H[14] - x356*H[6];
   const float x362 = -x355*H[15] - x356*H[7];
   const float x363 = -x355*H[16] - x356*H[8];
   const float x364 = -x355*H[11] - x356*H[3] + 1;
   const float x365 = x144 + x145 + x146;
   const float x366 = x136 + x138 + x140;
   const float x367 = x141 + x142 + x143;
   const float x368 = -x366*H[9] - x367*H[1];
   const float x369 = -x366*H[10] - x367*H[2];
   const float x370 = -x366*H[11] - x367*H[3];
   const float x371 = -x366*H[13] - x367*H[5];
   const float x372 = -x366*H[14] - x367*H[6];
   const float x373 = -x366*H[15] - x367*H[7];
   const float x374 = -x366*H[16] - x367*H[8];
   const float x375 = -x366*H[12] - x367*H[4] + 1;
   const float x376 = x156 + x157 + x158;
   const float x377 = x148 + x150 + x152;
   const float x378 = x153 + x154 + x155;
   const float x379 = -x377*H[9] - x378*H[1];
   const float x380 = -x377*H[10] - x378*H[2];
   const float x381 = -x377*H[11] - x378*H[3];
   const float x382 = -x377*H[12] - x378*H[4];
   const float x383 = -x377*H[14] - x378*H[6];
   const float x384 = -x377*H[15] - x378*H[7];
   const float x385 = -x377*H[16] - x378*H[8];
   const float x386 = -x377*H[13] - x378*H[5] + 1;
   const float x387 = x168 + x169 + x170;
   const float x388 = x160 + x162 + x164;
   const float x389 = x165 + x166 + x167;
   const float x390 = -x388*H[9] - x389*H[1];
   const float x391 = -x388*H[10] - x389*H[2];
   const float x392 = -x388*H[11] - x389*H[3];
   const float x393 = -x388*H[12] - x389*H[4];
   const float x394 = -x388*H[13] - x389*H[5];
   const float x395 = -x388*H[15] - x389*H[7];
   const float x396 = -x388*H[16] - x389*H[8];
   const float x397 = -x388*H[14] - x389*H[6] + 1;
   const float x398 = x180 + x181 + x182;
   const float x399 = x172 + x174 + x176;
   const float x400 = x177 + x178 + x179;
   const float x401 = -x399*H[9] - x400*H[1];
   const float x402 = -x399*H[10] - x400*H[2];
   const float x403 = -x399*H[11] - x400*H[3];
   const float x404 = -x399*H[12] - x400*H[4];
   const float x405 = -x399*H[13] - x400*H[5];
   const float x406 = -x399*H[14] - x400*H[6];
   const float x407 = -x399*H[16] - x400*H[8];
   const float x408 = -x399*H[15] - x400*H[7] + 1;
   const float x409 = x184 + x186;
   const float x410 = x187 + x188;
   const float x411 = -x409*H[9] - x410*H[1];
   const float x412 = -x409*H[10] - x410*H[2];
   const float x413 = -x409*H[11] - x410*H[3];
   const float x414 = -x409*H[12] - x410*H[4];
   const float x415 = -x409*H[13] - x410*H[5];
   const float x416 = -x409*H[14] - x410*H[6];
   const float x417 = -x409*H[15] - x410*H[7];
   const float x418 = -x409*H[16] - x410*H[8];
   const float x419 = x189 + x190;
   const float x420 = x192 + x194;
   const float x421 = x195 + x196;
   const float x422 = -x420*H[9] - x421*H[1];
   const float x423 = -x420*H[10] - x421*H[2];
   const float x424 = -x420*H[11] - x421*H[3];
   const float x425 = -x420*H[12] - x421*H[4];
   const float x426 = -x420*H[13] - x421*H[5];
   const float x427 = -x420*H[14] - x421*H[6];
   const float x428 = -x420*H[15] - x421*H[7];
   const float x429 = -x420*H[16] - x421*H[8];
   const float x430 = x197 + x198;
   const float x431 = x200 + x202;
   const float x432 = x203 + x204;
   const float x433 = -x431*H[9] - x432*H[1];
   const float x434 = -x431*H[10] - x432*H[2];
   const float x435 = -x431*H[11] - x432*H[3];
   const float x436 = -x431*H[12] - x432*H[4];
   const float x437 = -x431*H[13] - x432*H[5];
   const float x438 = -x431*H[14] - x432*H[6];
   const float x439 = -x431*H[15] - x432*H[7];
   const float x440 = -x431*H[16] - x432*H[8];
   const float x441 = x205 + x206;
   const float x442 = x209 + x211;
   const float x443 = x212 + x213;
   const float x444 = -x442*H[16] - x443*H[8];
   const float x445 = -x442*H[12] - x443*H[4];
   const float x446 = -x442*H[13] - x443*H[5];
   const float x447 = -x442*H[14] - x443*H[6];
   const float x448 = -x442*H[15] - x443*H[7];
   const float x449 = x214 + x215;
   const float x450 = x218 + x220;
   const float x451 = x221 + x222;
   const float x452 = -x450*H[16] - x451*H[8];
   const float x453 = -x450*H[12] - x451*H[4];
   const float x454 = -x450*H[13] - x451*H[5];
   const float x455 = -x450*H[14] - x451*H[6];
   const float x456 = -x450*H[15] - x451*H[7];
   const float x457 = x223 + x224;
   const float x458 = x227 + x229;
   const float x459 = x230 + x231;
   const float x460 = -x458*H[16] - x459*H[8];
   const float x461 = -x458*H[12] - x459*H[4];
   const float x462 = -x458*H[13] - x459*H[5];
   const float x463 = -x458*H[14] - x459*H[6];
   const float x464 = -x458*H[15] - x459*H[7];
   const float x465 = x232 + x233;
   const float x466 = x243 + x244 + x245;
   const float x467 = -x466*H[17] + 1;
   const float x468 = x235 + x237 + x239;
   const float x469 = x240 + x241 + x242;
   const float x470 = -x468*H[9] - x469*H[1];
   const float x471 = -x468*H[10] - x469*H[2];
   const float x472 = -x468*H[11] - x469*H[3];
   const float x473 = -x468*H[12] - x469*H[4];
   const float x474 = -x468*H[13] - x469*H[5];
   const float x475 = -x468*H[14] - x469*H[6];
   const float x476 = -x468*H[15] - x469*H[7];
   const float x477 = -x468*H[16] - x469*H[8];
   const float x478 = x248 + x250;
   const float x479 = x251 + x252;
   const float x480 = -x478*H[16] - x479*H[8];
   const float x481 = -x478*H[12] - x479*H[4];
   const float x482 = -x478*H[13] - x479*H[5];
   const float x483 = -x478*H[14] - x479*H[6];
   const float x484 = -x478*H[15] - x479*H[7];
   const float x485 = x253 + x254;
   const float x486 = x257 + x259;
   const float x487 = x260 + x261;
   const float x488 = -x486*H[16] - x487*H[8];
   const float x489 = -x486*H[12] - x487*H[4];
   const float x490 = -x486*H[13] - x487*H[5];
   const float x491 = -x486*H[14] - x487*H[6];
   const float x492 = -x486*H[15] - x487*H[7];
   const float x493 = x262 + x263;
   const float x494 = x266 + x268;
   const float x495 = x269 + x270;
   const float x496 = -x494*H[16] - x495*H[8];
   const float x497 = -x494*H[12] - x495*H[4];
   const float x498 = -x494*H[13] - x495*H[5];
   const float x499 = -x494*H[14] - x495*H[6];
   const float x500 = -x494*H[15] - x495*H[7];
   const float x501 = x271 + x272;
   const float x502 = x274 + x276;
   const float x503 = x277 + x278;
   const float x504 = -x502*H[12] - x503*H[4];
   const float x505 = -x502*H[13] - x503*H[5];
   const float x506 = -x502*H[14] - x503*H[6];
   const float x507 = -x502*H[15] - x503*H[7];
   const float x508 = -x502*H[16] - x503*H[8];
   const float x509 = x279 + x280;
   const float x510 = x282 + x284;
   const float x511 = x285 + x286;
   const float x512 = -x510*H[12] - x511*H[4];
   const float x513 = -x510*H[13] - x511*H[5];
   const float x514 = -x510*H[14] - x511*H[6];
   const float x515 = -x510*H[15] - x511*H[7];
   const float x516 = -x510*H[16] - x511*H[8];
   const float x517 = x287 + x288;
   const float x518 = x290 + x292;
   const float x519 = x293 + x294;
   const float x520 = -x518*H[12] - x519*H[4];
   const float x521 = -x518*H[13] - x519*H[5];
   const float x522 = -x518*H[14] - x519*H[6];
   const float x523 = -x518*H[15] - x519*H[7];
   const float x524 = -x518*H[16] - x519*H[8];
   const float x525 = x295 + x296;
   const float x526 = x306 + x307 + x308;
   const float x527 = -x526*H[18] + 1;
   const float x528 = x298 + x300 + x302;
   const float x529 = x303 + x304 + x305;
   const float x530 = -x528*H[9] - x529*H[1];
   const float x531 = -x528*H[10] - x529*H[2];
   const float x532 = -x528*H[11] - x529*H[3];
   const float x533 = -x528*H[12] - x529*H[4];
   const float x534 = -x528*H[13] - x529*H[5];
   const float x535 = -x528*H[14] - x529*H[6];
   const float x536 = -x528*H[15] - x529*H[7];
   const float x537 = -x528*H[16] - x529*H[8];
   const float x538 = x318 + x319 + x320;
   const float x539 = x310 + x312 + x314;
   const float x540 = x315 + x316 + x317;
   xnew[0] = x7*(-x76 - x79 - x84) + x85*(-x87 - x89 - x91) + x92*(-x94 - x96 - x98) + x[0];
   xnew[1] = x7*(-x100 - x102 - x104) + x85*(-x105 - x106 - x107) + x92*(-x108 - x109 - x110) + x[1];
   xnew[2] = x7*(-x112 - x114 - x116) + x85*(-x117 - x118 - x119) + x92*(-x120 - x121 - x122) + x[2];
   xnew[3] = x7*(-x124 - x126 - x128) + x85*(-x129 - x130 - x131) + x92*(-x132 - x133 - x134) + x[3];
   xnew[4] = x7*(-x136 - x138 - x140) + x85*(-x141 - x142 - x143) + x92*(-x144 - x145 - x146) + x[4];
   xnew[5] = x7*(-x148 - x150 - x152) + x85*(-x153 - x154 - x155) + x92*(-x156 - x157 - x158) + x[5];
   xnew[6] = x7*(-x160 - x162 - x164) + x85*(-x165 - x166 - x167) + x92*(-x168 - x169 - x170) + x[6];
   xnew[7] = x7*(-x172 - x174 - x176) + x85*(-x177 - x178 - x179) + x92*(-x180 - x181 - x182) + x[7];
   xnew[8] = x7*(-x184 - x186) + x85*(-x187 - x188) + x92*(-x189 - x190) + x[8];
   xnew[9] = x7*(-x192 - x194) + x85*(-x195 - x196) + x92*(-x197 - x198) + x[9];
   xnew[10] = x7*(-x200 - x202) + x85*(-x203 - x204) + x92*(-x205 - x206) + x[10];
   xnew[11] = x7*(-x209 - x211) + x85*(-x212 - x213) + x92*(-x214 - x215) + x[11];
   xnew[12] = x7*(-x218 - x220) + x85*(-x221 - x222) + x92*(-x223 - x224) + x[12];
   xnew[13] = x7*(-x227 - x229) + x85*(-x230 - x231) + x92*(-x232 - x233) + x[13];
   xnew[14] = x7*(-x235 - x237 - x239) + x85*(-x240 - x241 - x242) + x92*(-x243 - x244 - x245) + x[14];
   xnew[15] = x7*(-x248 - x250) + x85*(-x251 - x252) + x92*(-x253 - x254) + x[15];
   xnew[16] = x7*(-x257 - x259) + x85*(-x260 - x261) + x92*(-x262 - x263) + x[16];
   xnew[17] = x7*(-x266 - x268) + x85*(-x269 - x270) + x92*(-x271 - x272) + x[17];
   xnew[18] = x7*(-x274 - x276) + x85*(-x277 - x278) + x92*(-x279 - x280) + x[18];
   xnew[19] = x7*(-x282 - x284) + x85*(-x285 - x286) + x92*(-x287 - x288) + x[19];
   xnew[20] = x7*(-x290 - x292) + x85*(-x293 - x294) + x92*(-x295 - x296) + x[20];
   xnew[21] = x7*(-x298 - x300 - x302) + x85*(-x303 - x304 - x305) + x92*(-x306 - x307 - x308) + x[21];
   xnew[22] = x7*(-x310 - x312 - x314) + x85*(-x315 - x316 - x317) + x92*(-x318 - x319 - x320) + x[22];
   Pnew[0] = -x321*x81 - x321*x82 + x324*P[1] + x325*P[2] + x326*P[3] + x327*P[4] + x328*P[5] + x329*P[6] + x330*P[7] + x331*P[13] + P[0];
   Pnew[1] = -x11*x321 - x12*x321 + x324*P[14] + x325*P[15] + x326*P[16] + x327*P[17] + x328*P[18] + x329*P[19] + x330*P[20] + x331*P[26] + P[1];
   Pnew[2] = -x14*x321 - x15*x321 + x324*P[15] + x325*P[27] + x326*P[28] + x327*P[29] + x328*P[30] + x329*P[31] + x330*P[32] + x331*P[38] + P[2];
   Pnew[3] = -x17*x321 - x18*x321 + x324*P[16] + x325*P[28] + x326*P[39] + x327*P[40] + x328*P[41] + x329*P[42] + x330*P[43] + x331*P[49] + P[3];
   Pnew[4] = -x20*x321 - x21*x321 + x324*P[17] + x325*P[29] + x326*P[40] + x327*P[50] + x328*P[51] + x329*P[52] + x330*P[53] + x331*P[62] + P[4];
   Pnew[5] = -x23*x321 - x24*x321 + x324*P[18] + x325*P[30] + x326*P[41] + x327*P[51] + x328*P[63] + x329*P[64] + x330*P[65] + x331*P[74] + P[5];
   Pnew[6] = -x26*x321 - x27*x321 + x324*P[19] + x325*P[31] + x326*P[42] + x327*P[52] + x328*P[64] + x329*P[75] + x330*P[76] + x331*P[85] + P[6];
   Pnew[7] = -x29*x321 - x30*x321 + x324*P[20] + x325*P[32] + x326*P[43] + x327*P[53] + x328*P[65] + x329*P[76] + x330*P[86] + x331*P[95] + P[7];
   Pnew[8] = x324*P[21] + x325*P[33] + x326*P[44] + x327*P[54] + x328*P[66] + x329*P[77] + x330*P[87] + x331*P[101] + P[8];
   Pnew[9] = x324*P[22] + x325*P[34] + x326*P[45] + x327*P[55] + x328*P[67] + x329*P[78] + x330*P[88] + x331*P[105] + P[9];
   Pnew[10] = x324*P[23] + x325*P[35] + x326*P[46] + x327*P[56] + x328*P[68] + x329*P[79] + x330*P[89] + x331*P[110] + P[10];
   Pnew[11] = -x321*x36 - x321*x37 + x324*P[24] + x325*P[36] + x326*P[47] + x327*P[57] + x328*P[69] + x329*P[80] + x330*P[90] + x331*P[125] + P[11];
   Pnew[12] = -x321*x39 - x321*x40 + x324*P[25] + x325*P[37] + x326*P[48] + x327*P[61] + x328*P[73] + x329*P[84] + x330*P[94] + x331*P[142] + P[12];
   Pnew[13] = -x32*x321 - x321*x33 + x324*P[26] + x325*P[38] + x326*P[49] + x327*P[62] + x328*P[74] + x329*P[85] + x330*P[95] + x331*P[143] + P[13];
   Pnew[14] = -x11*x332 - x12*x332 + x335*P[15] + x336*P[16] + x337*P[17] + x338*P[18] + x339*P[19] + x340*P[20] + x341*P[26] + x342*P[14];
   Pnew[15] = -x14*x332 - x15*x332 + x335*P[27] + x336*P[28] + x337*P[29] + x338*P[30] + x339*P[31] + x340*P[32] + x341*P[38] + x342*P[15];
   Pnew[16] = -x17*x332 - x18*x332 + x335*P[28] + x336*P[39] + x337*P[40] + x338*P[41] + x339*P[42] + x340*P[43] + x341*P[49] + x342*P[16];
   Pnew[17] = -x20*x332 - x21*x332 + x335*P[29] + x336*P[40] + x337*P[50] + x338*P[51] + x339*P[52] + x340*P[53] + x341*P[62] + x342*P[17];
   Pnew[18] = -x23*x332 - x24*x332 + x335*P[30] + x336*P[41] + x337*P[51] + x338*P[63] + x339*P[64] + x340*P[65] + x341*P[74] + x342*P[18];
   Pnew[19] = -x26*x332 - x27*x332 + x335*P[31] + x336*P[42] + x337*P[52] + x338*P[64] + x339*P[75] + x340*P[76] + x341*P[85] + x342*P[19];
   Pnew[20] = -x29*x332 - x30*x332 + x335*P[32] + x336*P[43] + x337*P[53] + x338*P[65] + x339*P[76] + x340*P[86] + x341*P[95] + x342*P[20];
   Pnew[21] = x335*P[33] + x336*P[44] + x337*P[54] + x338*P[66] + x339*P[77] + x340*P[87] + x341*P[101] + x342*P[21];
   Pnew[22] = x335*P[34] + x336*P[45] + x337*P[55] + x338*P[67] + x339*P[78] + x340*P[88] + x341*P[105] + x342*P[22];
   Pnew[23] = x335*P[35] + x336*P[46] + x337*P[56] + x338*P[68] + x339*P[79] + x340*P[89] + x341*P[110] + x342*P[23];
   Pnew[24] = -x332*x36 - x332*x37 + x335*P[36] + x336*P[47] + x337*P[57] + x338*P[69] + x339*P[80] + x340*P[90] + x341*P[125] + x342*P[24];
   Pnew[25] = -x332*x39 - x332*x40 + x335*P[37] + x336*P[48] + x337*P[61] + x338*P[73] + x339*P[84] + x340*P[94] + x341*P[142] + x342*P[25];
   Pnew[26] = -x32*x332 - x33*x332 + x335*P[38] + x336*P[49] + x337*P[62] + x338*P[74] + x339*P[85] + x340*P[95] + x341*P[143] + x342*P[26];
   Pnew[27] = -x14*x343 - x15*x343 + x346*P[15] + x347*P[28] + x348*P[29] + x349*P[30] + x350*P[31] + x351*P[32] + x352*P[38] + x353*P[27];
   Pnew[28] = -x17*x343 - x18*x343 + x346*P[16] + x347*P[39] + x348*P[40] + x349*P[41] + x350*P[42] + x351*P[43] + x352*P[49] + x353*P[28];
   Pnew[29] = -x20*x343 - x21*x343 + x346*P[17] + x347*P[40] + x348*P[50] + x349*P[51] + x350*P[52] + x351*P[53] + x352*P[62] + x353*P[29];
   Pnew[30] = -x23*x343 - x24*x343 + x346*P[18] + x347*P[41] + x348*P[51] + x349*P[63] + x350*P[64] + x351*P[65] + x352*P[74] + x353*P[30];
   Pnew[31] = -x26*x343 - x27*x343 + x346*P[19] + x347*P[42] + x348*P[52] + x349*P[64] + x350*P[75] + x351*P[76] + x352*P[85] + x353*P[31];
   Pnew[32] = -x29*x343 - x30*x343 + x346*P[20] + x347*P[43] + x348*P[53] + x349*P[65] + x350*P[76] + x351*P[86] + x352*P[95] + x353*P[32];
   Pnew[33] = x346*P[21] + x347*P[44] + x348*P[54] + x349*P[66] + x350*P[77] + x351*P[87] + x352*P[101] + x353*P[33];
   Pnew[34] = x346*P[22] + x347*P[45] + x348*P[55] + x349*P[67] + x350*P[78] + x351*P[88] + x352*P[105] + x353*P[34];
   Pnew[35] = x346*P[23] + x347*P[46] + x348*P[56] + x349*P[68] + x350*P[79] + x351*P[89] + x352*P[110] + x353*P[35];
   Pnew[36] = -x343*x36 - x343*x37 + x346*P[24] + x347*P[47] + x348*P[57] + x349*P[69] + x350*P[80] + x351*P[90] + x352*P[125] + x353*P[36];
   Pnew[37] = -x343*x39 - x343*x40 + x346*P[25] + x347*P[48] + x348*P[61] + x349*P[73] + x350*P[84] + x351*P[94] + x352*P[142] + x353*P[37];
   Pnew[38] = -x32*x343 - x33*x343 + x346*P[26] + x347*P[49] + x348*P[62] + x349*P[74] + x350*P[85] + x351*P[95] + x352*P[143] + x353*P[38];
   Pnew[39] = -x17*x354 - x18*x354 + x357*P[16] + x358*P[28] + x359*P[40] + x360*P[41] + x361*P[42] + x362*P[43] + x363*P[49] + x364*P[39];
   Pnew[40] = -x20*x354 - x21*x354 + x357*P[17] + x358*P[29] + x359*P[50] + x360*P[51] + x361*P[52] + x362*P[53] + x363*P[62] + x364*P[40];
   Pnew[41] = -x23*x354 - x24*x354 + x357*P[18] + x358*P[30] + x359*P[51] + x360*P[63] + x361*P[64] + x362*P[65] + x363*P[74] + x364*P[41];
   Pnew[42] = -x26*x354 - x27*x354 + x357*P[19] + x358*P[31] + x359*P[52] + x360*P[64] + x361*P[75] + x362*P[76] + x363*P[85] + x364*P[42];
   Pnew[43] = -x29*x354 - x30*x354 + x357*P[20] + x358*P[32] + x359*P[53] + x360*P[65] + x361*P[76] + x362*P[86] + x363*P[95] + x364*P[43];
   Pnew[44] = x357*P[21] + x358*P[33] + x359*P[54] + x360*P[66] + x361*P[77] + x362*P[87] + x363*P[101] + x364*P[44];
   Pnew[45] = x357*P[22] + x358*P[34] + x359*P[55] + x360*P[67] + x361*P[78] + x362*P[88] + x363*P[105] + x364*P[45];
   Pnew[46] = x357*P[23] + x358*P[35] + x359*P[56] + x360*P[68] + x361*P[79] + x362*P[89] + x363*P[110] + x364*P[46];
   Pnew[47] = -x354*x36 - x354*x37 + x357*P[24] + x358*P[36] + x359*P[57] + x360*P[69] + x361*P[80] + x362*P[90] + x363*P[125] + x364*P[47];
   Pnew[48] = -x354*x39 - x354*x40 + x357*P[25] + x358*P[37] + x359*P[61] + x360*P[73] + x361*P[84] + x362*P[94] + x363*P[142] + x364*P[48];
   Pnew[49] = -x32*x354 - x33*x354 + x357*P[26] + x358*P[38] + x359*P[62] + x360*P[74] + x361*P[85] + x362*P[95] + x363*P[143] + x364*P[49];
   Pnew[50] = -x20*x365 - x21*x365 + x368*P[17] + x369*P[29] + x370*P[40] + x371*P[51] + x372*P[52] + x373*P[53] + x374*P[62] + x375*P[50];
   Pnew[51] = -x23*x365 - x24*x365 + x368*P[18] + x369*P[30] + x370*P[41] + x371*P[63] + x372*P[64] + x373*P[65] + x374*P[74] + x375*P[51];
   Pnew[52] = -x26*x365 - x27*x365 + x368*P[19] + x369*P[31] + x370*P[42] + x371*P[64] + x372*P[75] + x373*P[76] + x374*P[85] + x375*P[52];
   Pnew[53] = -x29*x365 - x30*x365 + x368*P[20] + x369*P[32] + x370*P[43] + x371*P[65] + x372*P[76] + x373*P[86] + x374*P[95] + x375*P[53];
   Pnew[54] = x368*P[21] + x369*P[33] + x370*P[44] + x371*P[66] + x372*P[77] + x373*P[87] + x374*P[101] + x375*P[54];
   Pnew[55] = x368*P[22] + x369*P[34] + x370*P[45] + x371*P[67] + x372*P[78] + x373*P[88] + x374*P[105] + x375*P[55];
   Pnew[56] = x368*P[23] + x369*P[35] + x370*P[46] + x371*P[68] + x372*P[79] + x373*P[89] + x374*P[110] + x375*P[56];
   Pnew[57] = -x36*x365 - x365*x37 + x368*P[24] + x369*P[36] + x370*P[47] + x371*P[69] + x372*P[80] + x373*P[90] + x374*P[125] + x375*P[57];
   Pnew[58] = x371*P[70] + x372*P[81] + x373*P[91] + x374*P[136] + x375*P[58];
   Pnew[59] = x371*P[71] + x372*P[82] + x373*P[92] + x374*P[138] + x375*P[59];
   Pnew[60] = x371*P[72] + x372*P[83] + x373*P[93] + x374*P[140] + x375*P[60];
   Pnew[61] = -x365*x39 - x365*x40 + x368*P[25] + x369*P[37] + x370*P[48] + x371*P[73] + x372*P[84] + x373*P[94] + x374*P[142] + x375*P[61];
   Pnew[62] = -x32*x365 - x33*x365 + x368*P[26] + x369*P[38] + x370*P[49] + x371*P[74] + x372*P[85] + x373*P[95] + x374*P[143] + x375*P[62];
   Pnew[63] = -x23*x376 - x24*x376 + x379*P[18] + x380*P[30] + x381*P[41] + x382*P[51] + x383*P[64] + x384*P[65] + x385*P[74] + x386*P[63];
   Pnew[64] = -x26*x376 - x27*x376 + x379*P[19] + x380*P[31] + x381*P[42] + x382*P[52] + x383*P[75] + x384*P[76] + x385*P[85] + x386*P[64];
   Pnew[65] = -x29*x376 - x30*x376 + x379*P[20] + x380*P[32] + x381*P[43] + x382*P[53] + x383*P[76] + x384*P[86] + x385*P[95] + x386*P[65];
   Pnew[66] = x379*P[21] + x380*P[33] + x381*P[44] + x382*P[54] + x383*P[77] + x384*P[87] + x385*P[101] + x386*P[66];
   Pnew[67] = x379*P[22] + x380*P[34] + x381*P[45] + x382*P[55] + x383*P[78] + x384*P[88] + x385*P[105] + x386*P[67];
   Pnew[68] = x379*P[23] + x380*P[35] + x381*P[46] + x382*P[56] + x383*P[79] + x384*P[89] + x385*P[110] + x386*P[68];
   Pnew[69] = -x36*x376 - x37*x376 + x379*P[24] + x380*P[36] + x381*P[47] + x382*P[57] + x383*P[80] + x384*P[90] + x385*P[125] + x386*P[69];
   Pnew[70] = x382*P[58] + x383*P[81] + x384*P[91] + x385*P[136] + x386*P[70];
   Pnew[71] = x382*P[59] + x383*P[82] + x384*P[92] + x385*P[138] + x386*P[71];
   Pnew[72] = x382*P[60] + x383*P[83] + x384*P[93] + x385*P[140] + x386*P[72];
   Pnew[73] = -x376*x39 - x376*x40 + x379*P[25] + x380*P[37] + x381*P[48] + x382*P[61] + x383*P[84] + x384*P[94] + x385*P[142] + x386*P[73];
   Pnew[74] = -x32*x376 - x33*x376 + x379*P[26] + x380*P[38] + x381*P[49] + x382*P[62] + x383*P[85] + x384*P[95] + x385*P[143] + x386*P[74];
   Pnew[75] = -x26*x387 - x27*x387 + x390*P[19] + x391*P[31] + x392*P[42] + x393*P[52] + x394*P[64] + x395*P[76] + x396*P[85] + x397*P[75];
   Pnew[76] = -x29*x387 - x30*x387 + x390*P[20] + x391*P[32] + x392*P[43] + x393*P[53] + x394*P[65] + x395*P[86] + x396*P[95] + x397*P[76];
   Pnew[77] = x390*P[21] + x391*P[33] + x392*P[44] + x393*P[54] + x394*P[66] + x395*P[87] + x396*P[101] + x397*P[77];
   Pnew[78] = x390*P[22] + x391*P[34] + x392*P[45] + x393*P[55] + x394*P[67] + x395*P[88] + x396*P[105] + x397*P[78];
   Pnew[79] = x390*P[23] + x391*P[35] + x392*P[46] + x393*P[56] + x394*P[68] + x395*P[89] + x396*P[110] + x397*P[79];
   Pnew[80] = -x36*x387 - x37*x387 + x390*P[24] + x391*P[36] + x392*P[47] + x393*P[57] + x394*P[69] + x395*P[90] + x396*P[125] + x397*P[80];
   Pnew[81] = x393*P[58] + x394*P[70] + x395*P[91] + x396*P[136] + x397*P[81];
   Pnew[82] = x393*P[59] + x394*P[71] + x395*P[92] + x396*P[138] + x397*P[82];
   Pnew[83] = x393*P[60] + x394*P[72] + x395*P[93] + x396*P[140] + x397*P[83];
   Pnew[84] = -x387*x39 - x387*x40 + x390*P[25] + x391*P[37] + x392*P[48] + x393*P[61] + x394*P[73] + x395*P[94] + x396*P[142] + x397*P[84];
   Pnew[85] = -x32*x387 - x33*x387 + x390*P[26] + x391*P[38] + x392*P[49] + x393*P[62] + x394*P[74] + x395*P[95] + x396*P[143] + x397*P[85];
   Pnew[86] = -x29*x398 - x30*x398 + x401*P[20] + x402*P[32] + x403*P[43] + x404*P[53] + x405*P[65] + x406*P[76] + x407*P[95] + x408*P[86];
   Pnew[87] = x401*P[21] + x402*P[33] + x403*P[44] + x404*P[54] + x405*P[66] + x406*P[77] + x407*P[101] + x408*P[87];
   Pnew[88] = x401*P[22] + x402*P[34] + x403*P[45] + x404*P[55] + x405*P[67] + x406*P[78] + x407*P[105] + x408*P[88];
   Pnew[89] = x401*P[23] + x402*P[35] + x403*P[46] + x404*P[56] + x405*P[68] + x406*P[79] + x407*P[110] + x408*P[89];
   Pnew[90] = -x36*x398 - x37*x398 + x401*P[24] + x402*P[36] + x403*P[47] + x404*P[57] + x405*P[69] + x406*P[80] + x407*P[125] + x408*P[90];
   Pnew[91] = x404*P[58] + x405*P[70] + x406*P[81] + x407*P[136] + x408*P[91];
   Pnew[92] = x404*P[59] + x405*P[71] + x406*P[82] + x407*P[138] + x408*P[92];
   Pnew[93] = x404*P[60] + x405*P[72] + x406*P[83] + x407*P[140] + x408*P[93];
   Pnew[94] = -x39*x398 - x398*x40 + x401*P[25] + x402*P[37] + x403*P[48] + x404*P[61] + x405*P[73] + x406*P[84] + x407*P[142] + x408*P[94];
   Pnew[95] = -x32*x398 - x33*x398 + x401*P[26] + x402*P[38] + x403*P[49] + x404*P[62] + x405*P[74] + x406*P[85] + x407*P[143] + x408*P[95];
   Pnew[96] = x411*P[21] + x412*P[33] + x413*P[44] + x414*P[54] + x415*P[66] + x416*P[77] + x417*P[87] + x418*P[101] + P[96];
   Pnew[97] = x418*P[114] + P[97];
   Pnew[98] = x418*P[118] + P[98];
   Pnew[99] = x418*P[128] + P[99];
   Pnew[100] = x414*P[58] + x415*P[70] + x416*P[81] + x417*P[91] + x418*P[136] + P[100];
   Pnew[101] = -x32*x419 - x33*x419 + x411*P[26] + x412*P[38] + x413*P[49] + x414*P[62] + x415*P[74] + x416*P[85] + x417*P[95] + x418*P[143] + P[101];
   Pnew[102] = x422*P[22] + x423*P[34] + x424*P[45] + x425*P[55] + x426*P[67] + x427*P[78] + x428*P[88] + x429*P[105] + P[102];
   Pnew[103] = x429*P[131] + P[103];
   Pnew[104] = x425*P[59] + x426*P[71] + x427*P[82] + x428*P[92] + x429*P[138] + P[104];
   Pnew[105] = -x32*x430 - x33*x430 + x422*P[26] + x423*P[38] + x424*P[49] + x425*P[62] + x426*P[74] + x427*P[85] + x428*P[95] + x429*P[143] + P[105];
   Pnew[106] = x433*P[23] + x434*P[35] + x435*P[46] + x436*P[56] + x437*P[68] + x438*P[79] + x439*P[89] + x440*P[110] + P[106];
   Pnew[107] = x440*P[122] + P[107];
   Pnew[108] = x440*P[134] + P[108];
   Pnew[109] = x436*P[60] + x437*P[72] + x438*P[83] + x439*P[93] + x440*P[140] + P[109];
   Pnew[110] = -x32*x441 - x33*x441 + x433*P[26] + x434*P[38] + x435*P[49] + x436*P[62] + x437*P[74] + x438*P[85] + x439*P[95] + x440*P[143] + P[110];
   Pnew[111] = x444*P[114] + P[111];
   Pnew[112] = x444*P[128] + P[112];
   Pnew[113] = x444*P[136] + x445*P[58] + x446*P[70] + x447*P[81] + x448*P[91] + P[113];
   Pnew[114] = -x32*x449 - x33*x449 + x444*P[143] + x445*P[62] + x446*P[74] + x447*P[85] + x448*P[95] + (-x442*H[9] - x443*H[1])*P[26] + (-x442*H[10] - x443*H[2])*P[38] + (-x442*H[11] - x443*H[3])*P[49] + P[114];
   Pnew[115] = x452*P[118] + P[115];
   Pnew[116] = x452*P[131] + P[116];
   Pnew[117] = x452*P[138] + x453*P[59] + x454*P[71] + x455*P[82] + x456*P[92] + P[117];
   Pnew[118] = -x32*x457 - x33*x457 + x452*P[143] + x453*P[62] + x454*P[74] + x455*P[85] + x456*P[95] + (-x450*H[9] - x451*H[1])*P[26] + (-x450*H[10] - x451*H[2])*P[38] + (-x450*H[11] - x451*H[3])*P[49] + P[118];
   Pnew[119] = x460*P[122] + P[119];
   Pnew[120] = x460*P[134] + P[120];
   Pnew[121] = x460*P[140] + x461*P[60] + x462*P[72] + x463*P[83] + x464*P[93] + P[121];
   Pnew[122] = -x32*x465 - x33*x465 + x460*P[143] + x461*P[62] + x462*P[74] + x463*P[85] + x464*P[95] + (-x458*H[9] - x459*H[1])*P[26] + (-x458*H[10] - x459*H[2])*P[38] + (-x458*H[11] - x459*H[3])*P[49] + P[122];
   Pnew[123] = -x37*x466 + x467*P[123] + x470*P[24] + x471*P[36] + x472*P[47] + x473*P[57] + x474*P[69] + x475*P[80] + x476*P[90] + x477*P[125];
   Pnew[124] = -x40*x466 + x467*P[124] + x470*P[25] + x471*P[37] + x472*P[48] + x473*P[61] + x474*P[73] + x475*P[84] + x476*P[94] + x477*P[142];
   Pnew[125] = -x33*x466 + x467*P[125] + x470*P[26] + x471*P[38] + x472*P[49] + x473*P[62] + x474*P[74] + x475*P[85] + x476*P[95] + x477*P[143];
   Pnew[126] = x480*P[128] + P[126];
   Pnew[127] = x480*P[136] + x481*P[58] + x482*P[70] + x483*P[81] + x484*P[91] + P[127];
   Pnew[128] = -x32*x485 - x33*x485 + x480*P[143] + x481*P[62] + x482*P[74] + x483*P[85] + x484*P[95] + (-x478*H[9] - x479*H[1])*P[26] + (-x478*H[10] - x479*H[2])*P[38] + (-x478*H[11] - x479*H[3])*P[49] + P[128];
   Pnew[129] = x488*P[131] + P[129];
   Pnew[130] = x488*P[138] + x489*P[59] + x490*P[71] + x491*P[82] + x492*P[92] + P[130];
   Pnew[131] = -x32*x493 - x33*x493 + x488*P[143] + x489*P[62] + x490*P[74] + x491*P[85] + x492*P[95] + (-x486*H[9] - x487*H[1])*P[26] + (-x486*H[10] - x487*H[2])*P[38] + (-x486*H[11] - x487*H[3])*P[49] + P[131];
   Pnew[132] = x496*P[134] + P[132];
   Pnew[133] = x496*P[140] + x497*P[60] + x498*P[72] + x499*P[83] + x500*P[93] + P[133];
   Pnew[134] = -x32*x501 - x33*x501 + x496*P[143] + x497*P[62] + x498*P[74] + x499*P[85] + x500*P[95] + (-x494*H[9] - x495*H[1])*P[26] + (-x494*H[10] - x495*H[2])*P[38] + (-x494*H[11] - x495*H[3])*P[49] + P[134];
   Pnew[135] = x504*P[58] + x505*P[70] + x506*P[81] + x507*P[91] + x508*P[136] + P[135];
   Pnew[136] = -x32*x509 - x33*x509 + x504*P[62] + x505*P[74] + x506*P[85] + x507*P[95] + x508*P[143] + (-x502*H[9] - x503*H[1])*P[26] + (-x502*H[10] - x503*H[2])*P[38] + (-x502*H[11] - x503*H[3])*P[49] + P[136];
   Pnew[137] = x512*P[59] + x513*P[71] + x514*P[82] + x515*P[92] + x516*P[138] + P[137];
   Pnew[138] = -x32*x517 - x33*x517 + x512*P[62] + x513*P[74] + x514*P[85] + x515*P[95] + x516*P[143] + (-x510*H[9] - x511*H[1])*P[26] + (-x510*H[10] - x511*H[2])*P[38] + (-x510*H[11] - x511*H[3])*P[49] + P[138];
   Pnew[139] = x520*P[60] + x521*P[72] + x522*P[83] + x523*P[93] + x524*P[140] + P[139];
   Pnew[140] = -x32*x525 - x33*x525 + x520*P[62] + x521*P[74] + x522*P[85] + x523*P[95] + x524*P[143] + (-x518*H[9] - x519*H[1])*P[26] + (-x518*H[10] - x519*H[2])*P[38] + (-x518*H[11] - x519*H[3])*P[49] + P[140];
   Pnew[141] = -x39*x526 + x527*P[141] + x530*P[25] + x531*P[37] + x532*P[48] + x533*P[61] + x534*P[73] + x535*P[84] + x536*P[94] + x537*P[142];
   Pnew[142] = -x32*x526 + x527*P[142] + x530*P[26] + x531*P[38] + x532*P[49] + x533*P[62] + x534*P[74] + x535*P[85] + x536*P[95] + x537*P[143];
   Pnew[143] = -x32*x538 - x33*x538 + (-x539*H[9] - x540*H[1])*P[26] + (-x539*H[10] - x540*H[2])*P[38] + (-x539*H[11] - x540*H[3])*P[49] + (-x539*H[12] - x540*H[4])*P[62] + (-x539*H[13] - x540*H[5])*P[74] + (-x539*H[14] - x540*H[6])*P[85] + (-x539*H[15] - x540*H[7])*P[95] + (-x539*H[16] - x540*H[8] + 1)*P[143];
}

void covariance_init(float *Pnew) {
   float P0[NUMX] = {1.0e7f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0f,
                     1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                     1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                     1.0e-9f, 1.0e-9f, 1.0e-9f, 1.0e-9f,
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
   Pnew[13] = 0;
   Pnew[14] = P0[1];
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
   Pnew[25] = 0;
   Pnew[26] = 0;
   Pnew[27] = P0[2];
   Pnew[28] = 0;
   Pnew[29] = 0;
   Pnew[30] = 0;
   Pnew[31] = 0;
   Pnew[32] = 0;
   Pnew[33] = 0;
   Pnew[34] = 0;
   Pnew[35] = 0;
   Pnew[36] = 0;
   Pnew[37] = 0;
   Pnew[38] = 0;
   Pnew[39] = P0[3];
   Pnew[40] = 0;
   Pnew[41] = 0;
   Pnew[42] = 0;
   Pnew[43] = 0;
   Pnew[44] = 0;
   Pnew[45] = 0;
   Pnew[46] = 0;
   Pnew[47] = 0;
   Pnew[48] = 0;
   Pnew[49] = 0;
   Pnew[50] = P0[4];
   Pnew[51] = 0;
   Pnew[52] = 0;
   Pnew[53] = 0;
   Pnew[54] = 0;
   Pnew[55] = 0;
   Pnew[56] = 0;
   Pnew[57] = 0;
   Pnew[58] = 0;
   Pnew[59] = 0;
   Pnew[60] = 0;
   Pnew[61] = 0;
   Pnew[62] = 0;
   Pnew[63] = P0[5];
   Pnew[64] = 0;
   Pnew[65] = 0;
   Pnew[66] = 0;
   Pnew[67] = 0;
   Pnew[68] = 0;
   Pnew[69] = 0;
   Pnew[70] = 0;
   Pnew[71] = 0;
   Pnew[72] = 0;
   Pnew[73] = 0;
   Pnew[74] = 0;
   Pnew[75] = P0[6];
   Pnew[76] = 0;
   Pnew[77] = 0;
   Pnew[78] = 0;
   Pnew[79] = 0;
   Pnew[80] = 0;
   Pnew[81] = 0;
   Pnew[82] = 0;
   Pnew[83] = 0;
   Pnew[84] = 0;
   Pnew[85] = 0;
   Pnew[86] = P0[7];
   Pnew[87] = 0;
   Pnew[88] = 0;
   Pnew[89] = 0;
   Pnew[90] = 0;
   Pnew[91] = 0;
   Pnew[92] = 0;
   Pnew[93] = 0;
   Pnew[94] = 0;
   Pnew[95] = 0;
   Pnew[96] = P0[8];
   Pnew[97] = 0;
   Pnew[98] = 0;
   Pnew[99] = 0;
   Pnew[100] = 0;
   Pnew[101] = 0;
   Pnew[102] = P0[9];
   Pnew[103] = 0;
   Pnew[104] = 0;
   Pnew[105] = 0;
   Pnew[106] = P0[10];
   Pnew[107] = 0;
   Pnew[108] = 0;
   Pnew[109] = 0;
   Pnew[110] = 0;
   Pnew[111] = P0[11];
   Pnew[112] = 0;
   Pnew[113] = 0;
   Pnew[114] = 0;
   Pnew[115] = P0[12];
   Pnew[116] = 0;
   Pnew[117] = 0;
   Pnew[118] = 0;
   Pnew[119] = P0[13];
   Pnew[120] = 0;
   Pnew[121] = 0;
   Pnew[122] = 0;
   Pnew[123] = P0[14];
   Pnew[124] = 0;
   Pnew[125] = 0;
   Pnew[126] = P0[15];
   Pnew[127] = 0;
   Pnew[128] = 0;
   Pnew[129] = P0[16];
   Pnew[130] = 0;
   Pnew[131] = 0;
   Pnew[132] = P0[17];
   Pnew[133] = 0;
   Pnew[134] = 0;
   Pnew[135] = P0[18];
   Pnew[136] = 0;
   Pnew[137] = P0[19];
   Pnew[138] = 0;
   Pnew[139] = P0[20];
   Pnew[140] = 0;
   Pnew[141] = P0[21];
   Pnew[142] = 0;
   Pnew[143] = P0[22];
}

/**
 * @}
 * @}
 */
