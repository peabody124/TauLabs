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
	qcins_state->x[22] = 1.0f; // Initial mu to 1

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
	Q[18] = 1e-5f;                      // Thrust gain
	Q[19] = Q[20] = Q[21] = 1e-5f;      // Output bias states
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
	qcins_state->init.beta_t = 2.0f; // thrust to weight ratio 2:1
	qcins_state->init.mu = 1.0f;

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

bool qcins_get_thrust(uintptr_t qcins_handle, float thrust[1])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	thrust[0] = qcins_state->x[18];
	return true;
}

bool qcins_get_output_bias(uintptr_t qcins_handle, float out_bias[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	out_bias[0] = qcins_state->x[19];
	out_bias[1] = qcins_state->x[20];
	out_bias[2] = qcins_state->x[21];
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
   const float x2 = param[0]*x[14]*x[21];
   const float x3 = 2*x2;
   const float x4 = x[4]*x[7];
   const float x5 = x[5]*x[6];
   const float x6 = 2*x4 - 2*x5;
   const float x7 = x[4]*x[5];
   const float x8 = x[6]*x[7];
   const float x9 = x7 + x8;
   const float x10 = 2*x[3];
   const float x11 = (x[6]*x[6]);
   const float x12 = (x[5]*x[5]);
   const float x13 = -x12;
   const float x14 = (x[4]*x[4]);
   const float x15 = (x[7]*x[7]);
   const float x16 = x14 - x15;
   const float x17 = x11 + x13 + x16;
   const float x18 = (x10*x9 + x17*x[2] - x6*x[1])*x[22];
   const float x19 = -x11;
   const float x20 = x12 + x16 + x19;
   const float x21 = 2*x4 + 2*x5;
   const float x22 = x0 - x1;
   const float x23 = (-x10*x22 + x20*x[1] + x21*x[2])*x[22];
   const float x24 = (1.0F/2.0F)*x[8];
   const float x25 = (1.0F/2.0F)*x[6];
   const float x26 = (1.0F/2.0F)*x[7];
   const float x27 = (1.0F/2.0F)*x[4];
   const float x28 = (1.0F/2.0F)*x[5];
   const float x29 = Ts/param[5];
   xnew[0] = Ts*x[3] + x[0];
   xnew[1] = Ts*(x18*x6 - x20*x23 - x3*(x0 + x1)) + x[1];
   xnew[2] = Ts*(-x17*x18 - x21*x23 + x3*(x7 - x8)) + x[2];
   xnew[3] = Ts*(-2*x18*x9 - x2*(x13 + x14 + x15 + x19) + 2*x22*x23 + param[0]) + x[3];
   xnew[4] = Ts*(-x24*x[5] - x25*x[9] - x26*x[10]) + x[4];
   xnew[5] = Ts*(x24*x[4] + x25*x[10] - x26*x[9]) + x[5];
   xnew[6] = Ts*(x24*x[7] + x27*x[9] - x28*x[10]) + x[6];
   xnew[7] = Ts*(-x24*x[6] + x27*x[10] + x28*x[9]) + x[7];
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
   xnew[19] = x[19];
   xnew[20] = x[20];
   xnew[21] = x[21];
   xnew[22] = x[22];
}

void linearize_FH(const float * restrict x, const float Ts, const float * restrict param, float * restrict F, float * restrict H)
{
   const float x0 = 2*x[4];
   const float x1 = x0*x[7];
   const float x2 = 2*x[5];
   const float x3 = x2*x[6];
   const float x4 = -x1 + x3;
   const float x5 = (x[5]*x[5]);
   const float x6 = (x[6]*x[6]);
   const float x7 = -x6;
   const float x8 = (x[4]*x[4]);
   const float x9 = (x[7]*x[7]);
   const float x10 = x8 - x9;
   const float x11 = x10 + x5 + x7;
   const float x12 = -x5;
   const float x13 = x10 + x12 + x6;
   const float x14 = x4*x[22];
   const float x15 = x1 + x3;
   const float x16 = x11*x[22];
   const float x17 = Ts*(-x13*x14 - x15*x16);
   const float x18 = x0*x[5];
   const float x19 = 2*x[7];
   const float x20 = x19*x[6];
   const float x21 = x18 + x20;
   const float x22 = x0*x[6];
   const float x23 = x2*x[7];
   const float x24 = -x22 + x23;
   const float x25 = Ts*(-x14*x21 - x16*x24);
   const float x26 = x0*x[2] - x19*x[1] + x2*x[3];
   const float x27 = x26*x[22];
   const float x28 = x0*x[1];
   const float x29 = x19*x[2];
   const float x30 = 2*x[6];
   const float x31 = x30*x[3];
   const float x32 = x28 + x29 - x31;
   const float x33 = x24*x[3];
   const float x34 = x15*x[2];
   const float x35 = x11*x[1];
   const float x36 = x33 + x34 + x35;
   const float x37 = x36*x[22];
   const float x38 = x0*x37;
   const float x39 = x21*x[3];
   const float x40 = x4*x[1];
   const float x41 = x13*x[2];
   const float x42 = x39 + x40 + x41;
   const float x43 = x42*x[22];
   const float x44 = x19*x43;
   const float x45 = param[0]*x[21];
   const float x46 = x45*x[14];
   const float x47 = x30*x46;
   const float x48 = -x38 + x44 - x47;
   const float x49 = x30*x[1];
   const float x50 = x2*x[2];
   const float x51 = x0*x[3];
   const float x52 = x49 - x50 + x51;
   const float x53 = x19*x[3] + x2*x[1] + x30*x[2];
   const float x54 = -x19*x46 - x2*x37 - x30*x43;
   const float x55 = x53*x[22];
   const float x56 = -x49 + x50 - x51;
   const float x57 = x2*x43;
   const float x58 = x30*x37;
   const float x59 = x0*x46;
   const float x60 = -x57 + x58 - x59;
   const float x61 = x2*x46;
   const float x62 = -x28 - x29 + x31;
   const float x63 = x0*x43;
   const float x64 = x19*x37;
   const float x65 = Ts*(x22 + x23);
   const float x66 = param[0]*x[14];
   const float x67 = x15*x[22];
   const float x68 = x13*x[22];
   const float x69 = Ts*(-x21*x68 - x24*x67);
   const float x70 = x61 - x63 - x64;
   const float x71 = Ts*(-x18 + x20);
   const float x72 = x24*x[22];
   const float x73 = x21*x[22];
   const float x74 = Ts*(x12 + x7 + x8 + x9);
   const float x75 = (1.0F/2.0F)*Ts;
   const float x76 = x75*x[8];
   const float x77 = -x76;
   const float x78 = x75*x[9];
   const float x79 = -x78;
   const float x80 = x75*x[10];
   const float x81 = -x80;
   const float x82 = x75*x[5];
   const float x83 = -x82;
   const float x84 = x75*x[6];
   const float x85 = -x84;
   const float x86 = x75*x[7];
   const float x87 = -x86;
   const float x88 = x75*x[4];
   const float x89 = -Ts/param[5];
   const float x90 = -x55;
   const float x91 = -x27;
   const float x92 = x0*param[6];
   const float x93 = x19*param[7];
   const float x94 = x30*param[8];
   const float x95 = x92 + x93 - x94;
   const float x96 = x19*param[8] + x2*param[6] + x30*param[7];
   const float x97 = x30*param[6];
   const float x98 = x2*param[7];
   const float x99 = x0*param[8];
   const float x100 = x19*param[6];
   const float x101 = x0*param[7];
   const float x102 = x2*param[8];
   const float x103 = -x100 + x101 + x102;
   const float x104 = x97 - x98 + x99;
   F[0] = Ts;
   F[1] = Ts*(-(x11*x11)*x[22] - (x4*x4)*x[22]);
   F[2] = x17;
   F[3] = x25;
   F[4] = Ts*(-x16*x32 - x27*x4 + x48);
   F[5] = Ts*(-x14*x52 - x16*x53 + x54);
   F[6] = Ts*(-x16*x56 - x4*x55 + x60);
   F[7] = Ts*(-x14*x62 - x16*x26 - x61 + x63 + x64);
   F[8] = -x45*x65;
   F[9] = -x65*x66;
   F[10] = Ts*(-x11*x36 - x4*x42);
   F[11] = x17;
   F[12] = Ts*(-(x13*x13)*x[22] - (x15*x15)*x[22]);
   F[13] = x69;
   F[14] = Ts*(-x13*x27 - x32*x67 + x70);
   F[15] = Ts*(-x52*x68 - x53*x67 + x57 - x58 + x59);
   F[16] = Ts*(-x13*x55 + x54 - x56*x67);
   F[17] = Ts*(-x26*x67 + x48 - x62*x68);
   F[18] = -x45*x71;
   F[19] = -x66*x71;
   F[20] = Ts*(-x13*x42 - x15*x36);
   F[21] = x25;
   F[22] = x69;
   F[23] = Ts*(-(x21*x21)*x[22] - (x24*x24)*x[22]);
   F[24] = Ts*(-x21*x27 - x32*x72 + x60);
   F[25] = Ts*(-x52*x73 - x53*x72 + x70);
   F[26] = Ts*(-x21*x55 + x38 - x44 + x47 - x56*x72);
   F[27] = Ts*(-x26*x72 + x54 - x62*x73);
   F[28] = -x45*x74;
   F[29] = -x66*x74;
   F[30] = Ts*(-x21*x42 - x24*x36);
   F[31] = x77;
   F[32] = x79;
   F[33] = x81;
   F[34] = x83;
   F[35] = x85;
   F[36] = x87;
   F[37] = x76;
   F[38] = x80;
   F[39] = x79;
   F[40] = x88;
   F[41] = x87;
   F[42] = x84;
   F[43] = x78;
   F[44] = x81;
   F[45] = x76;
   F[46] = x86;
   F[47] = x88;
   F[48] = x83;
   F[49] = x80;
   F[50] = x78;
   F[51] = x77;
   F[52] = x85;
   F[53] = x82;
   F[54] = x88;
   F[55] = Ts*param[1];
   F[56] = Ts*param[2];
   F[57] = Ts*(param[3] - param[4]);
   F[58] = -Ts*param[4];
   F[59] = x89;
   F[60] = x89;
   F[61] = x89;
   F[62] = x89;
   F[63] = x89;
   F[64] = x89;
   F[65] = x89;
   H[0] = 1;
   H[1] = -x16;
   H[2] = -x67;
   H[3] = -x72;
   H[4] = -x32*x[22];
   H[5] = x90;
   H[6] = -x56*x[22];
   H[7] = x91;
   H[8] = -x33 - x34 - x35;
   H[9] = -x14;
   H[10] = -x68;
   H[11] = -x73;
   H[12] = x91;
   H[13] = -x52*x[22];
   H[14] = x90;
   H[15] = -x62*x[22];
   H[16] = -x39 - x40 - x41;
   H[17] = -x45;
   H[18] = -x66;
   H[19] = 1;
   H[20] = 1;
   H[21] = 1;
   H[22] = 1;
   H[23] = 1;
   H[24] = 1;
   H[25] = x95;
   H[26] = x96;
   H[27] = -x97 + x98 - x99;
   H[28] = x103;
   H[29] = x103;
   H[30] = x104;
   H[31] = x96;
   H[32] = -x92 - x93 + x94;
   H[33] = x104;
   H[34] = x100 - x101 - x102;
   H[35] = x95;
   H[36] = x96;
}

// TODO: rework F definition so that Ts is already multiplied into all the terms
// to reduce the memory requirement (currently duplicates F)
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
   const float x0 = 2*x[3];
   const float x1 = x[4]*x[7];
   const float x2 = x[5]*x[6];
   const float x3 = (x[6]*x[6]);
   const float x4 = (x[5]*x[5]);
   const float x5 = (x[4]*x[4]) - (x[7]*x[7]);
   const float x6 = -(x0*(x[4]*x[5] + x[6]*x[7]) - 2*(x1 - x2)*x[1] + (x3 - x4 + x5)*x[2])*x[22] - accel[1];
   const float x7 = H[1]*P[24] + H[2]*P[36] + H[3]*P[47] + H[4]*P[57] + H[5]*P[69] + H[6]*P[80] + H[7]*P[90] + H[8]*P[125];
   const float x8 = H[1]*P[25] + H[2]*P[37] + H[3]*P[48] + H[4]*P[61] + H[5]*P[73] + H[6]*P[84] + H[7]*P[94] + H[8]*P[142];
   const float x9 = x7*H[17] + x8*H[18];
   const float x10 = H[17]*P[24];
   const float x11 = H[18]*P[25];
   const float x12 = x10 + x11;
   const float x13 = H[17]*P[36];
   const float x14 = H[18]*P[37];
   const float x15 = x13 + x14;
   const float x16 = H[17]*P[47];
   const float x17 = H[18]*P[48];
   const float x18 = x16 + x17;
   const float x19 = H[17]*P[57];
   const float x20 = H[18]*P[61];
   const float x21 = x19 + x20;
   const float x22 = H[17]*P[69];
   const float x23 = H[18]*P[73];
   const float x24 = x22 + x23;
   const float x25 = H[17]*P[80];
   const float x26 = H[18]*P[84];
   const float x27 = x25 + x26;
   const float x28 = H[17]*P[90];
   const float x29 = H[18]*P[94];
   const float x30 = x28 + x29;
   const float x31 = H[17]*P[125];
   const float x32 = H[18]*P[142];
   const float x33 = x31 + x32;
   const float x34 = x12*H[9] + x15*H[10] + x18*H[11] + x21*H[12] + x24*H[13] + x27*H[14] + x30*H[15] + x33*H[16];
   const float x35 = H[17]*P[123];
   const float x36 = H[18]*P[124];
   const float x37 = x35 + x36;
   const float x38 = H[17]*P[124];
   const float x39 = H[18]*P[141];
   const float x40 = x38 + x39;
   const float x41 = x37*H[17] + x40*H[18] + R[3];
   const float x42 = H[1]*P[14] + H[2]*P[15] + H[3]*P[16] + H[4]*P[17] + H[5]*P[18] + H[6]*P[19] + H[7]*P[20] + H[8]*P[26];
   const float x43 = H[1]*P[15] + H[2]*P[27] + H[3]*P[28] + H[4]*P[29] + H[5]*P[30] + H[6]*P[31] + H[7]*P[32] + H[8]*P[38];
   const float x44 = H[1]*P[16] + H[2]*P[28] + H[3]*P[39] + H[4]*P[40] + H[5]*P[41] + H[6]*P[42] + H[7]*P[43] + H[8]*P[49];
   const float x45 = H[1]*P[17] + H[2]*P[29] + H[3]*P[40] + H[4]*P[50] + H[5]*P[51] + H[6]*P[52] + H[7]*P[53] + H[8]*P[62];
   const float x46 = H[1]*P[18] + H[2]*P[30] + H[3]*P[41] + H[4]*P[51] + H[5]*P[63] + H[6]*P[64] + H[7]*P[65] + H[8]*P[74];
   const float x47 = H[1]*P[19] + H[2]*P[31] + H[3]*P[42] + H[4]*P[52] + H[5]*P[64] + H[6]*P[75] + H[7]*P[76] + H[8]*P[85];
   const float x48 = H[1]*P[20] + H[2]*P[32] + H[3]*P[43] + H[4]*P[53] + H[5]*P[65] + H[6]*P[76] + H[7]*P[86] + H[8]*P[95];
   const float x49 = H[1]*P[26] + H[2]*P[38] + H[3]*P[49] + H[4]*P[62] + H[5]*P[74] + H[6]*P[85] + H[7]*P[95] + H[8]*P[143];
   const float x50 = x42*H[9] + x43*H[10] + x44*H[11] + x45*H[12] + x46*H[13] + x47*H[14] + x48*H[15] + x49*H[16];
   const float x51 = x34*x9 - x41*x50;
   const float x52 = H[9]*P[14] + H[10]*P[15] + H[11]*P[16] + H[12]*P[17] + H[13]*P[18] + H[14]*P[19] + H[15]*P[20] + H[16]*P[26];
   const float x53 = H[9]*P[15] + H[10]*P[27] + H[11]*P[28] + H[12]*P[29] + H[13]*P[30] + H[14]*P[31] + H[15]*P[32] + H[16]*P[38];
   const float x54 = H[9]*P[16] + H[10]*P[28] + H[11]*P[39] + H[12]*P[40] + H[13]*P[41] + H[14]*P[42] + H[15]*P[43] + H[16]*P[49];
   const float x55 = H[9]*P[17] + H[10]*P[29] + H[11]*P[40] + H[12]*P[50] + H[13]*P[51] + H[14]*P[52] + H[15]*P[53] + H[16]*P[62];
   const float x56 = H[9]*P[18] + H[10]*P[30] + H[11]*P[41] + H[12]*P[51] + H[13]*P[63] + H[14]*P[64] + H[15]*P[65] + H[16]*P[74];
   const float x57 = H[9]*P[19] + H[10]*P[31] + H[11]*P[42] + H[12]*P[52] + H[13]*P[64] + H[14]*P[75] + H[15]*P[76] + H[16]*P[85];
   const float x58 = H[9]*P[20] + H[10]*P[32] + H[11]*P[43] + H[12]*P[53] + H[13]*P[65] + H[14]*P[76] + H[15]*P[86] + H[16]*P[95];
   const float x59 = H[9]*P[26] + H[10]*P[38] + H[11]*P[49] + H[12]*P[62] + H[13]*P[74] + H[14]*P[85] + H[15]*P[95] + H[16]*P[143];
   const float x60 = x52*H[1] + x53*H[2] + x54*H[3] + x55*H[4] + x56*H[5] + x57*H[6] + x58*H[7] + x59*H[8];
   const float x61 = x34*x60;
   const float x62 = H[9]*P[24] + H[10]*P[36] + H[11]*P[47] + H[12]*P[57] + H[13]*P[69] + H[14]*P[80] + H[15]*P[90] + H[16]*P[125];
   const float x63 = H[9]*P[25] + H[10]*P[37] + H[11]*P[48] + H[12]*P[61] + H[13]*P[73] + H[14]*P[84] + H[15]*P[94] + H[16]*P[142];
   const float x64 = x62*H[17] + x63*H[18];
   const float x65 = x12*H[1] + x15*H[2] + x18*H[3] + x21*H[4] + x24*H[5] + x27*H[6] + x30*H[7] + x33*H[8];
   const float x66 = x64*x65;
   const float x67 = x52*H[9] + x53*H[10] + x54*H[11] + x55*H[12] + x56*H[13] + x57*H[14] + x58*H[15] + x59*H[16] + R[2];
   const float x68 = x65*x67;
   const float x69 = x42*H[1] + x43*H[2] + x44*H[3] + x45*H[4] + x46*H[5] + x47*H[6] + x48*H[7] + x49*H[8] + R[1];
   const float x70 = x34*x64;
   const float x71 = x41*x60;
   const float x72 = x41*x67;
   const float x73 = 1.0F/(x50*x66 - x50*x71 + x61*x9 - x68*x9 - x69*x70 + x69*x72);
   const float x74 = x73*(H[1]*P[1] + H[2]*P[2] + H[3]*P[3] + H[4]*P[4] + H[5]*P[5] + H[6]*P[6] + H[7]*P[7] + H[8]*P[13]);
   const float x75 = x51*x74;
   const float x76 = x41*x69 - x65*x9;
   const float x77 = x73*(H[9]*P[1] + H[10]*P[2] + H[11]*P[3] + H[12]*P[4] + H[13]*P[5] + H[14]*P[6] + H[15]*P[7] + H[16]*P[13]);
   const float x78 = x76*x77;
   const float x79 = -x34*x69 + x50*x65;
   const float x80 = H[17]*P[11];
   const float x81 = H[18]*P[12];
   const float x82 = x73*(x80 + x81);
   const float x83 = x79*x82;
   const float x84 = -(-x0*(x[4]*x[6] - x[5]*x[7]) + 2*(x1 + x2)*x[2] + (-x3 + x4 + x5)*x[1])*x[22] - accel[0];
   const float x85 = x66 - x71;
   const float x86 = x77*x85;
   const float x87 = -x70 + x72;
   const float x88 = x74*x87;
   const float x89 = x61 - x68;
   const float x90 = x82*x89;
   const float x91 = -accel[2] - param[0]*x[14]*x[21];
   const float x92 = x60*x9 - x64*x69;
   const float x93 = x77*x92;
   const float x94 = x50*x64 - x67*x9;
   const float x95 = x74*x94;
   const float x96 = -x50*x60 + x67*x69;
   const float x97 = x82*x96;
   const float x98 = x42*x73;
   const float x99 = x51*x98;
   const float x100 = x52*x73;
   const float x101 = x100*x76;
   const float x102 = x12*x73;
   const float x103 = x102*x79;
   const float x104 = x100*x85;
   const float x105 = x87*x98;
   const float x106 = x102*x89;
   const float x107 = x100*x92;
   const float x108 = x94*x98;
   const float x109 = x102*x96;
   const float x110 = x43*x73;
   const float x111 = x110*x51;
   const float x112 = x53*x73;
   const float x113 = x112*x76;
   const float x114 = x15*x73;
   const float x115 = x114*x79;
   const float x116 = x112*x85;
   const float x117 = x110*x87;
   const float x118 = x114*x89;
   const float x119 = x112*x92;
   const float x120 = x110*x94;
   const float x121 = x114*x96;
   const float x122 = x44*x73;
   const float x123 = x122*x51;
   const float x124 = x54*x73;
   const float x125 = x124*x76;
   const float x126 = x18*x73;
   const float x127 = x126*x79;
   const float x128 = x124*x85;
   const float x129 = x122*x87;
   const float x130 = x126*x89;
   const float x131 = x124*x92;
   const float x132 = x122*x94;
   const float x133 = x126*x96;
   const float x134 = x45*x73;
   const float x135 = x134*x51;
   const float x136 = x55*x73;
   const float x137 = x136*x76;
   const float x138 = x21*x73;
   const float x139 = x138*x79;
   const float x140 = x136*x85;
   const float x141 = x134*x87;
   const float x142 = x138*x89;
   const float x143 = x136*x92;
   const float x144 = x134*x94;
   const float x145 = x138*x96;
   const float x146 = x46*x73;
   const float x147 = x146*x51;
   const float x148 = x56*x73;
   const float x149 = x148*x76;
   const float x150 = x24*x73;
   const float x151 = x150*x79;
   const float x152 = x148*x85;
   const float x153 = x146*x87;
   const float x154 = x150*x89;
   const float x155 = x148*x92;
   const float x156 = x146*x94;
   const float x157 = x150*x96;
   const float x158 = x47*x73;
   const float x159 = x158*x51;
   const float x160 = x57*x73;
   const float x161 = x160*x76;
   const float x162 = x27*x73;
   const float x163 = x162*x79;
   const float x164 = x160*x85;
   const float x165 = x158*x87;
   const float x166 = x162*x89;
   const float x167 = x160*x92;
   const float x168 = x158*x94;
   const float x169 = x162*x96;
   const float x170 = x48*x73;
   const float x171 = x170*x51;
   const float x172 = x58*x73;
   const float x173 = x172*x76;
   const float x174 = x30*x73;
   const float x175 = x174*x79;
   const float x176 = x172*x85;
   const float x177 = x170*x87;
   const float x178 = x174*x89;
   const float x179 = x172*x92;
   const float x180 = x170*x94;
   const float x181 = x174*x96;
   const float x182 = x73*(H[1]*P[21] + H[2]*P[33] + H[3]*P[44] + H[4]*P[54] + H[5]*P[66] + H[6]*P[77] + H[7]*P[87] + H[8]*P[101]);
   const float x183 = x182*x51;
   const float x184 = x73*(H[9]*P[21] + H[10]*P[33] + H[11]*P[44] + H[12]*P[54] + H[13]*P[66] + H[14]*P[77] + H[15]*P[87] + H[16]*P[101]);
   const float x185 = x184*x76;
   const float x186 = x184*x85;
   const float x187 = x182*x87;
   const float x188 = x184*x92;
   const float x189 = x182*x94;
   const float x190 = x73*(H[1]*P[22] + H[2]*P[34] + H[3]*P[45] + H[4]*P[55] + H[5]*P[67] + H[6]*P[78] + H[7]*P[88] + H[8]*P[105]);
   const float x191 = x190*x51;
   const float x192 = x73*(H[9]*P[22] + H[10]*P[34] + H[11]*P[45] + H[12]*P[55] + H[13]*P[67] + H[14]*P[78] + H[15]*P[88] + H[16]*P[105]);
   const float x193 = x192*x76;
   const float x194 = x192*x85;
   const float x195 = x190*x87;
   const float x196 = x192*x92;
   const float x197 = x190*x94;
   const float x198 = x73*(H[1]*P[23] + H[2]*P[35] + H[3]*P[46] + H[4]*P[56] + H[5]*P[68] + H[6]*P[79] + H[7]*P[89] + H[8]*P[110]);
   const float x199 = x198*x51;
   const float x200 = x73*(H[9]*P[23] + H[10]*P[35] + H[11]*P[46] + H[12]*P[56] + H[13]*P[68] + H[14]*P[79] + H[15]*P[89] + H[16]*P[110]);
   const float x201 = x200*x76;
   const float x202 = x200*x85;
   const float x203 = x198*x87;
   const float x204 = x200*x92;
   const float x205 = x198*x94;
   const float x206 = x73*P[114];
   const float x207 = x206*H[8];
   const float x208 = x207*x51;
   const float x209 = x206*H[16];
   const float x210 = x209*x76;
   const float x211 = x209*x85;
   const float x212 = x207*x87;
   const float x213 = x209*x92;
   const float x214 = x207*x94;
   const float x215 = x73*P[118];
   const float x216 = x215*H[8];
   const float x217 = x216*x51;
   const float x218 = x215*H[16];
   const float x219 = x218*x76;
   const float x220 = x218*x85;
   const float x221 = x216*x87;
   const float x222 = x218*x92;
   const float x223 = x216*x94;
   const float x224 = x73*P[122];
   const float x225 = x224*H[8];
   const float x226 = x225*x51;
   const float x227 = x224*H[16];
   const float x228 = x227*x76;
   const float x229 = x227*x85;
   const float x230 = x225*x87;
   const float x231 = x227*x92;
   const float x232 = x225*x94;
   const float x233 = x7*x73;
   const float x234 = x233*x51;
   const float x235 = x62*x73;
   const float x236 = x235*x76;
   const float x237 = x37*x73;
   const float x238 = x237*x79;
   const float x239 = x235*x85;
   const float x240 = x233*x87;
   const float x241 = x237*x89;
   const float x242 = x235*x92;
   const float x243 = x233*x94;
   const float x244 = x237*x96;
   const float x245 = x73*P[128];
   const float x246 = x245*H[8];
   const float x247 = x246*x51;
   const float x248 = x245*H[16];
   const float x249 = x248*x76;
   const float x250 = x248*x85;
   const float x251 = x246*x87;
   const float x252 = x248*x92;
   const float x253 = x246*x94;
   const float x254 = x73*P[131];
   const float x255 = x254*H[8];
   const float x256 = x255*x51;
   const float x257 = x254*H[16];
   const float x258 = x257*x76;
   const float x259 = x257*x85;
   const float x260 = x255*x87;
   const float x261 = x257*x92;
   const float x262 = x255*x94;
   const float x263 = x73*P[134];
   const float x264 = x263*H[8];
   const float x265 = x264*x51;
   const float x266 = x263*H[16];
   const float x267 = x266*x76;
   const float x268 = x266*x85;
   const float x269 = x264*x87;
   const float x270 = x266*x92;
   const float x271 = x264*x94;
   const float x272 = x73*(H[4]*P[58] + H[5]*P[70] + H[6]*P[81] + H[7]*P[91] + H[8]*P[136]);
   const float x273 = x272*x51;
   const float x274 = x73*(H[12]*P[58] + H[13]*P[70] + H[14]*P[81] + H[15]*P[91] + H[16]*P[136]);
   const float x275 = x274*x76;
   const float x276 = x274*x85;
   const float x277 = x272*x87;
   const float x278 = x274*x92;
   const float x279 = x272*x94;
   const float x280 = x73*(H[4]*P[59] + H[5]*P[71] + H[6]*P[82] + H[7]*P[92] + H[8]*P[138]);
   const float x281 = x280*x51;
   const float x282 = x73*(H[12]*P[59] + H[13]*P[71] + H[14]*P[82] + H[15]*P[92] + H[16]*P[138]);
   const float x283 = x282*x76;
   const float x284 = x282*x85;
   const float x285 = x280*x87;
   const float x286 = x282*x92;
   const float x287 = x280*x94;
   const float x288 = x73*(H[4]*P[60] + H[5]*P[72] + H[6]*P[83] + H[7]*P[93] + H[8]*P[140]);
   const float x289 = x288*x51;
   const float x290 = x73*(H[12]*P[60] + H[13]*P[72] + H[14]*P[83] + H[15]*P[93] + H[16]*P[140]);
   const float x291 = x290*x76;
   const float x292 = x290*x85;
   const float x293 = x288*x87;
   const float x294 = x290*x92;
   const float x295 = x288*x94;
   const float x296 = x73*x8;
   const float x297 = x296*x51;
   const float x298 = x63*x73;
   const float x299 = x298*x76;
   const float x300 = x40*x73;
   const float x301 = x300*x79;
   const float x302 = x298*x85;
   const float x303 = x296*x87;
   const float x304 = x300*x89;
   const float x305 = x298*x92;
   const float x306 = x296*x94;
   const float x307 = x300*x96;
   const float x308 = x49*x73;
   const float x309 = x308*x51;
   const float x310 = x59*x73;
   const float x311 = x310*x76;
   const float x312 = x33*x73;
   const float x313 = x312*x79;
   const float x314 = x310*x85;
   const float x315 = x308*x87;
   const float x316 = x312*x89;
   const float x317 = x310*x92;
   const float x318 = x308*x94;
   const float x319 = x312*x96;
   const float x320 = x93 + x95 + x97;
   const float x321 = x75 + x78 + x83;
   const float x322 = x86 + x88 + x90;
   const float x323 = -x321*H[9] - x322*H[1];
   const float x324 = -x321*H[10] - x322*H[2];
   const float x325 = -x321*H[11] - x322*H[3];
   const float x326 = -x321*H[12] - x322*H[4];
   const float x327 = -x321*H[13] - x322*H[5];
   const float x328 = -x321*H[14] - x322*H[6];
   const float x329 = -x321*H[15] - x322*H[7];
   const float x330 = -x321*H[16] - x322*H[8];
   const float x331 = x107 + x108 + x109;
   const float x332 = x101 + x103 + x99;
   const float x333 = x104 + x105 + x106;
   const float x334 = -x332*H[10] - x333*H[2];
   const float x335 = -x332*H[11] - x333*H[3];
   const float x336 = -x332*H[12] - x333*H[4];
   const float x337 = -x332*H[13] - x333*H[5];
   const float x338 = -x332*H[14] - x333*H[6];
   const float x339 = -x332*H[15] - x333*H[7];
   const float x340 = -x332*H[16] - x333*H[8];
   const float x341 = -x332*H[9] - x333*H[1] + 1;
   const float x342 = x119 + x120 + x121;
   const float x343 = x111 + x113 + x115;
   const float x344 = x116 + x117 + x118;
   const float x345 = -x343*H[9] - x344*H[1];
   const float x346 = -x343*H[11] - x344*H[3];
   const float x347 = -x343*H[12] - x344*H[4];
   const float x348 = -x343*H[13] - x344*H[5];
   const float x349 = -x343*H[14] - x344*H[6];
   const float x350 = -x343*H[15] - x344*H[7];
   const float x351 = -x343*H[16] - x344*H[8];
   const float x352 = -x343*H[10] - x344*H[2] + 1;
   const float x353 = x131 + x132 + x133;
   const float x354 = x123 + x125 + x127;
   const float x355 = x128 + x129 + x130;
   const float x356 = -x354*H[9] - x355*H[1];
   const float x357 = -x354*H[10] - x355*H[2];
   const float x358 = -x354*H[12] - x355*H[4];
   const float x359 = -x354*H[13] - x355*H[5];
   const float x360 = -x354*H[14] - x355*H[6];
   const float x361 = -x354*H[15] - x355*H[7];
   const float x362 = -x354*H[16] - x355*H[8];
   const float x363 = -x354*H[11] - x355*H[3] + 1;
   const float x364 = x143 + x144 + x145;
   const float x365 = x135 + x137 + x139;
   const float x366 = x140 + x141 + x142;
   const float x367 = -x365*H[9] - x366*H[1];
   const float x368 = -x365*H[10] - x366*H[2];
   const float x369 = -x365*H[11] - x366*H[3];
   const float x370 = -x365*H[13] - x366*H[5];
   const float x371 = -x365*H[14] - x366*H[6];
   const float x372 = -x365*H[15] - x366*H[7];
   const float x373 = -x365*H[16] - x366*H[8];
   const float x374 = -x365*H[12] - x366*H[4] + 1;
   const float x375 = x155 + x156 + x157;
   const float x376 = x147 + x149 + x151;
   const float x377 = x152 + x153 + x154;
   const float x378 = -x376*H[9] - x377*H[1];
   const float x379 = -x376*H[10] - x377*H[2];
   const float x380 = -x376*H[11] - x377*H[3];
   const float x381 = -x376*H[12] - x377*H[4];
   const float x382 = -x376*H[14] - x377*H[6];
   const float x383 = -x376*H[15] - x377*H[7];
   const float x384 = -x376*H[16] - x377*H[8];
   const float x385 = -x376*H[13] - x377*H[5] + 1;
   const float x386 = x167 + x168 + x169;
   const float x387 = x159 + x161 + x163;
   const float x388 = x164 + x165 + x166;
   const float x389 = -x387*H[9] - x388*H[1];
   const float x390 = -x387*H[10] - x388*H[2];
   const float x391 = -x387*H[11] - x388*H[3];
   const float x392 = -x387*H[12] - x388*H[4];
   const float x393 = -x387*H[13] - x388*H[5];
   const float x394 = -x387*H[15] - x388*H[7];
   const float x395 = -x387*H[16] - x388*H[8];
   const float x396 = -x387*H[14] - x388*H[6] + 1;
   const float x397 = x179 + x180 + x181;
   const float x398 = x171 + x173 + x175;
   const float x399 = x176 + x177 + x178;
   const float x400 = -x398*H[9] - x399*H[1];
   const float x401 = -x398*H[10] - x399*H[2];
   const float x402 = -x398*H[11] - x399*H[3];
   const float x403 = -x398*H[12] - x399*H[4];
   const float x404 = -x398*H[13] - x399*H[5];
   const float x405 = -x398*H[14] - x399*H[6];
   const float x406 = -x398*H[16] - x399*H[8];
   const float x407 = -x398*H[15] - x399*H[7] + 1;
   const float x408 = x183 + x185;
   const float x409 = x186 + x187;
   const float x410 = -x408*H[9] - x409*H[1];
   const float x411 = -x408*H[10] - x409*H[2];
   const float x412 = -x408*H[11] - x409*H[3];
   const float x413 = -x408*H[12] - x409*H[4];
   const float x414 = -x408*H[13] - x409*H[5];
   const float x415 = -x408*H[14] - x409*H[6];
   const float x416 = -x408*H[15] - x409*H[7];
   const float x417 = -x408*H[16] - x409*H[8];
   const float x418 = x188 + x189;
   const float x419 = x191 + x193;
   const float x420 = x194 + x195;
   const float x421 = -x419*H[9] - x420*H[1];
   const float x422 = -x419*H[10] - x420*H[2];
   const float x423 = -x419*H[11] - x420*H[3];
   const float x424 = -x419*H[12] - x420*H[4];
   const float x425 = -x419*H[13] - x420*H[5];
   const float x426 = -x419*H[14] - x420*H[6];
   const float x427 = -x419*H[15] - x420*H[7];
   const float x428 = -x419*H[16] - x420*H[8];
   const float x429 = x196 + x197;
   const float x430 = x199 + x201;
   const float x431 = x202 + x203;
   const float x432 = -x430*H[9] - x431*H[1];
   const float x433 = -x430*H[10] - x431*H[2];
   const float x434 = -x430*H[11] - x431*H[3];
   const float x435 = -x430*H[12] - x431*H[4];
   const float x436 = -x430*H[13] - x431*H[5];
   const float x437 = -x430*H[14] - x431*H[6];
   const float x438 = -x430*H[15] - x431*H[7];
   const float x439 = -x430*H[16] - x431*H[8];
   const float x440 = x204 + x205;
   const float x441 = x208 + x210;
   const float x442 = x211 + x212;
   const float x443 = -x441*H[16] - x442*H[8];
   const float x444 = -x441*H[12] - x442*H[4];
   const float x445 = -x441*H[13] - x442*H[5];
   const float x446 = -x441*H[14] - x442*H[6];
   const float x447 = -x441*H[15] - x442*H[7];
   const float x448 = x213 + x214;
   const float x449 = x217 + x219;
   const float x450 = x220 + x221;
   const float x451 = -x449*H[16] - x450*H[8];
   const float x452 = -x449*H[12] - x450*H[4];
   const float x453 = -x449*H[13] - x450*H[5];
   const float x454 = -x449*H[14] - x450*H[6];
   const float x455 = -x449*H[15] - x450*H[7];
   const float x456 = x222 + x223;
   const float x457 = x226 + x228;
   const float x458 = x229 + x230;
   const float x459 = -x457*H[16] - x458*H[8];
   const float x460 = -x457*H[12] - x458*H[4];
   const float x461 = -x457*H[13] - x458*H[5];
   const float x462 = -x457*H[14] - x458*H[6];
   const float x463 = -x457*H[15] - x458*H[7];
   const float x464 = x231 + x232;
   const float x465 = x242 + x243 + x244;
   const float x466 = -x465*H[17] + 1;
   const float x467 = x234 + x236 + x238;
   const float x468 = x239 + x240 + x241;
   const float x469 = -x467*H[9] - x468*H[1];
   const float x470 = -x467*H[10] - x468*H[2];
   const float x471 = -x467*H[11] - x468*H[3];
   const float x472 = -x467*H[12] - x468*H[4];
   const float x473 = -x467*H[13] - x468*H[5];
   const float x474 = -x467*H[14] - x468*H[6];
   const float x475 = -x467*H[15] - x468*H[7];
   const float x476 = -x467*H[16] - x468*H[8];
   const float x477 = x247 + x249;
   const float x478 = x250 + x251;
   const float x479 = -x477*H[16] - x478*H[8];
   const float x480 = -x477*H[12] - x478*H[4];
   const float x481 = -x477*H[13] - x478*H[5];
   const float x482 = -x477*H[14] - x478*H[6];
   const float x483 = -x477*H[15] - x478*H[7];
   const float x484 = x252 + x253;
   const float x485 = x256 + x258;
   const float x486 = x259 + x260;
   const float x487 = -x485*H[16] - x486*H[8];
   const float x488 = -x485*H[12] - x486*H[4];
   const float x489 = -x485*H[13] - x486*H[5];
   const float x490 = -x485*H[14] - x486*H[6];
   const float x491 = -x485*H[15] - x486*H[7];
   const float x492 = x261 + x262;
   const float x493 = x265 + x267;
   const float x494 = x268 + x269;
   const float x495 = -x493*H[16] - x494*H[8];
   const float x496 = -x493*H[12] - x494*H[4];
   const float x497 = -x493*H[13] - x494*H[5];
   const float x498 = -x493*H[14] - x494*H[6];
   const float x499 = -x493*H[15] - x494*H[7];
   const float x500 = x270 + x271;
   const float x501 = x273 + x275;
   const float x502 = x276 + x277;
   const float x503 = -x501*H[12] - x502*H[4];
   const float x504 = -x501*H[13] - x502*H[5];
   const float x505 = -x501*H[14] - x502*H[6];
   const float x506 = -x501*H[15] - x502*H[7];
   const float x507 = -x501*H[16] - x502*H[8];
   const float x508 = x278 + x279;
   const float x509 = x281 + x283;
   const float x510 = x284 + x285;
   const float x511 = -x509*H[12] - x510*H[4];
   const float x512 = -x509*H[13] - x510*H[5];
   const float x513 = -x509*H[14] - x510*H[6];
   const float x514 = -x509*H[15] - x510*H[7];
   const float x515 = -x509*H[16] - x510*H[8];
   const float x516 = x286 + x287;
   const float x517 = x289 + x291;
   const float x518 = x292 + x293;
   const float x519 = -x517*H[12] - x518*H[4];
   const float x520 = -x517*H[13] - x518*H[5];
   const float x521 = -x517*H[14] - x518*H[6];
   const float x522 = -x517*H[15] - x518*H[7];
   const float x523 = -x517*H[16] - x518*H[8];
   const float x524 = x294 + x295;
   const float x525 = x305 + x306 + x307;
   const float x526 = -x525*H[18] + 1;
   const float x527 = x297 + x299 + x301;
   const float x528 = x302 + x303 + x304;
   const float x529 = -x527*H[9] - x528*H[1];
   const float x530 = -x527*H[10] - x528*H[2];
   const float x531 = -x527*H[11] - x528*H[3];
   const float x532 = -x527*H[12] - x528*H[4];
   const float x533 = -x527*H[13] - x528*H[5];
   const float x534 = -x527*H[14] - x528*H[6];
   const float x535 = -x527*H[15] - x528*H[7];
   const float x536 = -x527*H[16] - x528*H[8];
   const float x537 = x317 + x318 + x319;
   const float x538 = x309 + x311 + x313;
   const float x539 = x314 + x315 + x316;
   xnew[0] = x6*(-x75 - x78 - x83) + x84*(-x86 - x88 - x90) + x91*(-x93 - x95 - x97) + x[0];
   xnew[1] = x6*(-x101 - x103 - x99) + x84*(-x104 - x105 - x106) + x91*(-x107 - x108 - x109) + x[1];
   xnew[2] = x6*(-x111 - x113 - x115) + x84*(-x116 - x117 - x118) + x91*(-x119 - x120 - x121) + x[2];
   xnew[3] = x6*(-x123 - x125 - x127) + x84*(-x128 - x129 - x130) + x91*(-x131 - x132 - x133) + x[3];
   xnew[4] = x6*(-x135 - x137 - x139) + x84*(-x140 - x141 - x142) + x91*(-x143 - x144 - x145) + x[4];
   xnew[5] = x6*(-x147 - x149 - x151) + x84*(-x152 - x153 - x154) + x91*(-x155 - x156 - x157) + x[5];
   xnew[6] = x6*(-x159 - x161 - x163) + x84*(-x164 - x165 - x166) + x91*(-x167 - x168 - x169) + x[6];
   xnew[7] = x6*(-x171 - x173 - x175) + x84*(-x176 - x177 - x178) + x91*(-x179 - x180 - x181) + x[7];
   xnew[8] = x6*(-x183 - x185) + x84*(-x186 - x187) + x91*(-x188 - x189) + x[8];
   xnew[9] = x6*(-x191 - x193) + x84*(-x194 - x195) + x91*(-x196 - x197) + x[9];
   xnew[10] = x6*(-x199 - x201) + x84*(-x202 - x203) + x91*(-x204 - x205) + x[10];
   xnew[11] = x6*(-x208 - x210) + x84*(-x211 - x212) + x91*(-x213 - x214) + x[11];
   xnew[12] = x6*(-x217 - x219) + x84*(-x220 - x221) + x91*(-x222 - x223) + x[12];
   xnew[13] = x6*(-x226 - x228) + x84*(-x229 - x230) + x91*(-x231 - x232) + x[13];
   xnew[14] = x6*(-x234 - x236 - x238) + x84*(-x239 - x240 - x241) + x91*(-x242 - x243 - x244) + x[14];
   xnew[15] = x6*(-x247 - x249) + x84*(-x250 - x251) + x91*(-x252 - x253) + x[15];
   xnew[16] = x6*(-x256 - x258) + x84*(-x259 - x260) + x91*(-x261 - x262) + x[16];
   xnew[17] = x6*(-x265 - x267) + x84*(-x268 - x269) + x91*(-x270 - x271) + x[17];
   xnew[18] = x6*(-x273 - x275) + x84*(-x276 - x277) + x91*(-x278 - x279) + x[18];
   xnew[19] = x6*(-x281 - x283) + x84*(-x284 - x285) + x91*(-x286 - x287) + x[19];
   xnew[20] = x6*(-x289 - x291) + x84*(-x292 - x293) + x91*(-x294 - x295) + x[20];
   xnew[21] = x6*(-x297 - x299 - x301) + x84*(-x302 - x303 - x304) + x91*(-x305 - x306 - x307) + x[21];
   xnew[22] = x6*(-x309 - x311 - x313) + x84*(-x314 - x315 - x316) + x91*(-x317 - x318 - x319) + x[22];
   Pnew[0] = -x320*x80 - x320*x81 + x323*P[1] + x324*P[2] + x325*P[3] + x326*P[4] + x327*P[5] + x328*P[6] + x329*P[7] + x330*P[13] + P[0];
   Pnew[1] = -x10*x320 - x11*x320 + x323*P[14] + x324*P[15] + x325*P[16] + x326*P[17] + x327*P[18] + x328*P[19] + x329*P[20] + x330*P[26] + P[1];
   Pnew[2] = -x13*x320 - x14*x320 + x323*P[15] + x324*P[27] + x325*P[28] + x326*P[29] + x327*P[30] + x328*P[31] + x329*P[32] + x330*P[38] + P[2];
   Pnew[3] = -x16*x320 - x17*x320 + x323*P[16] + x324*P[28] + x325*P[39] + x326*P[40] + x327*P[41] + x328*P[42] + x329*P[43] + x330*P[49] + P[3];
   Pnew[4] = -x19*x320 - x20*x320 + x323*P[17] + x324*P[29] + x325*P[40] + x326*P[50] + x327*P[51] + x328*P[52] + x329*P[53] + x330*P[62] + P[4];
   Pnew[5] = -x22*x320 - x23*x320 + x323*P[18] + x324*P[30] + x325*P[41] + x326*P[51] + x327*P[63] + x328*P[64] + x329*P[65] + x330*P[74] + P[5];
   Pnew[6] = -x25*x320 - x26*x320 + x323*P[19] + x324*P[31] + x325*P[42] + x326*P[52] + x327*P[64] + x328*P[75] + x329*P[76] + x330*P[85] + P[6];
   Pnew[7] = -x28*x320 - x29*x320 + x323*P[20] + x324*P[32] + x325*P[43] + x326*P[53] + x327*P[65] + x328*P[76] + x329*P[86] + x330*P[95] + P[7];
   Pnew[8] = x323*P[21] + x324*P[33] + x325*P[44] + x326*P[54] + x327*P[66] + x328*P[77] + x329*P[87] + x330*P[101] + P[8];
   Pnew[9] = x323*P[22] + x324*P[34] + x325*P[45] + x326*P[55] + x327*P[67] + x328*P[78] + x329*P[88] + x330*P[105] + P[9];
   Pnew[10] = x323*P[23] + x324*P[35] + x325*P[46] + x326*P[56] + x327*P[68] + x328*P[79] + x329*P[89] + x330*P[110] + P[10];
   Pnew[11] = -x320*x35 - x320*x36 + x323*P[24] + x324*P[36] + x325*P[47] + x326*P[57] + x327*P[69] + x328*P[80] + x329*P[90] + x330*P[125] + P[11];
   Pnew[12] = -x320*x38 - x320*x39 + x323*P[25] + x324*P[37] + x325*P[48] + x326*P[61] + x327*P[73] + x328*P[84] + x329*P[94] + x330*P[142] + P[12];
   Pnew[13] = -x31*x320 - x32*x320 + x323*P[26] + x324*P[38] + x325*P[49] + x326*P[62] + x327*P[74] + x328*P[85] + x329*P[95] + x330*P[143] + P[13];
   Pnew[14] = -x10*x331 - x11*x331 + x334*P[15] + x335*P[16] + x336*P[17] + x337*P[18] + x338*P[19] + x339*P[20] + x340*P[26] + x341*P[14];
   Pnew[15] = -x13*x331 - x14*x331 + x334*P[27] + x335*P[28] + x336*P[29] + x337*P[30] + x338*P[31] + x339*P[32] + x340*P[38] + x341*P[15];
   Pnew[16] = -x16*x331 - x17*x331 + x334*P[28] + x335*P[39] + x336*P[40] + x337*P[41] + x338*P[42] + x339*P[43] + x340*P[49] + x341*P[16];
   Pnew[17] = -x19*x331 - x20*x331 + x334*P[29] + x335*P[40] + x336*P[50] + x337*P[51] + x338*P[52] + x339*P[53] + x340*P[62] + x341*P[17];
   Pnew[18] = -x22*x331 - x23*x331 + x334*P[30] + x335*P[41] + x336*P[51] + x337*P[63] + x338*P[64] + x339*P[65] + x340*P[74] + x341*P[18];
   Pnew[19] = -x25*x331 - x26*x331 + x334*P[31] + x335*P[42] + x336*P[52] + x337*P[64] + x338*P[75] + x339*P[76] + x340*P[85] + x341*P[19];
   Pnew[20] = -x28*x331 - x29*x331 + x334*P[32] + x335*P[43] + x336*P[53] + x337*P[65] + x338*P[76] + x339*P[86] + x340*P[95] + x341*P[20];
   Pnew[21] = x334*P[33] + x335*P[44] + x336*P[54] + x337*P[66] + x338*P[77] + x339*P[87] + x340*P[101] + x341*P[21];
   Pnew[22] = x334*P[34] + x335*P[45] + x336*P[55] + x337*P[67] + x338*P[78] + x339*P[88] + x340*P[105] + x341*P[22];
   Pnew[23] = x334*P[35] + x335*P[46] + x336*P[56] + x337*P[68] + x338*P[79] + x339*P[89] + x340*P[110] + x341*P[23];
   Pnew[24] = -x331*x35 - x331*x36 + x334*P[36] + x335*P[47] + x336*P[57] + x337*P[69] + x338*P[80] + x339*P[90] + x340*P[125] + x341*P[24];
   Pnew[25] = -x331*x38 - x331*x39 + x334*P[37] + x335*P[48] + x336*P[61] + x337*P[73] + x338*P[84] + x339*P[94] + x340*P[142] + x341*P[25];
   Pnew[26] = -x31*x331 - x32*x331 + x334*P[38] + x335*P[49] + x336*P[62] + x337*P[74] + x338*P[85] + x339*P[95] + x340*P[143] + x341*P[26];
   Pnew[27] = -x13*x342 - x14*x342 + x345*P[15] + x346*P[28] + x347*P[29] + x348*P[30] + x349*P[31] + x350*P[32] + x351*P[38] + x352*P[27];
   Pnew[28] = -x16*x342 - x17*x342 + x345*P[16] + x346*P[39] + x347*P[40] + x348*P[41] + x349*P[42] + x350*P[43] + x351*P[49] + x352*P[28];
   Pnew[29] = -x19*x342 - x20*x342 + x345*P[17] + x346*P[40] + x347*P[50] + x348*P[51] + x349*P[52] + x350*P[53] + x351*P[62] + x352*P[29];
   Pnew[30] = -x22*x342 - x23*x342 + x345*P[18] + x346*P[41] + x347*P[51] + x348*P[63] + x349*P[64] + x350*P[65] + x351*P[74] + x352*P[30];
   Pnew[31] = -x25*x342 - x26*x342 + x345*P[19] + x346*P[42] + x347*P[52] + x348*P[64] + x349*P[75] + x350*P[76] + x351*P[85] + x352*P[31];
   Pnew[32] = -x28*x342 - x29*x342 + x345*P[20] + x346*P[43] + x347*P[53] + x348*P[65] + x349*P[76] + x350*P[86] + x351*P[95] + x352*P[32];
   Pnew[33] = x345*P[21] + x346*P[44] + x347*P[54] + x348*P[66] + x349*P[77] + x350*P[87] + x351*P[101] + x352*P[33];
   Pnew[34] = x345*P[22] + x346*P[45] + x347*P[55] + x348*P[67] + x349*P[78] + x350*P[88] + x351*P[105] + x352*P[34];
   Pnew[35] = x345*P[23] + x346*P[46] + x347*P[56] + x348*P[68] + x349*P[79] + x350*P[89] + x351*P[110] + x352*P[35];
   Pnew[36] = -x342*x35 - x342*x36 + x345*P[24] + x346*P[47] + x347*P[57] + x348*P[69] + x349*P[80] + x350*P[90] + x351*P[125] + x352*P[36];
   Pnew[37] = -x342*x38 - x342*x39 + x345*P[25] + x346*P[48] + x347*P[61] + x348*P[73] + x349*P[84] + x350*P[94] + x351*P[142] + x352*P[37];
   Pnew[38] = -x31*x342 - x32*x342 + x345*P[26] + x346*P[49] + x347*P[62] + x348*P[74] + x349*P[85] + x350*P[95] + x351*P[143] + x352*P[38];
   Pnew[39] = -x16*x353 - x17*x353 + x356*P[16] + x357*P[28] + x358*P[40] + x359*P[41] + x360*P[42] + x361*P[43] + x362*P[49] + x363*P[39];
   Pnew[40] = -x19*x353 - x20*x353 + x356*P[17] + x357*P[29] + x358*P[50] + x359*P[51] + x360*P[52] + x361*P[53] + x362*P[62] + x363*P[40];
   Pnew[41] = -x22*x353 - x23*x353 + x356*P[18] + x357*P[30] + x358*P[51] + x359*P[63] + x360*P[64] + x361*P[65] + x362*P[74] + x363*P[41];
   Pnew[42] = -x25*x353 - x26*x353 + x356*P[19] + x357*P[31] + x358*P[52] + x359*P[64] + x360*P[75] + x361*P[76] + x362*P[85] + x363*P[42];
   Pnew[43] = -x28*x353 - x29*x353 + x356*P[20] + x357*P[32] + x358*P[53] + x359*P[65] + x360*P[76] + x361*P[86] + x362*P[95] + x363*P[43];
   Pnew[44] = x356*P[21] + x357*P[33] + x358*P[54] + x359*P[66] + x360*P[77] + x361*P[87] + x362*P[101] + x363*P[44];
   Pnew[45] = x356*P[22] + x357*P[34] + x358*P[55] + x359*P[67] + x360*P[78] + x361*P[88] + x362*P[105] + x363*P[45];
   Pnew[46] = x356*P[23] + x357*P[35] + x358*P[56] + x359*P[68] + x360*P[79] + x361*P[89] + x362*P[110] + x363*P[46];
   Pnew[47] = -x35*x353 - x353*x36 + x356*P[24] + x357*P[36] + x358*P[57] + x359*P[69] + x360*P[80] + x361*P[90] + x362*P[125] + x363*P[47];
   Pnew[48] = -x353*x38 - x353*x39 + x356*P[25] + x357*P[37] + x358*P[61] + x359*P[73] + x360*P[84] + x361*P[94] + x362*P[142] + x363*P[48];
   Pnew[49] = -x31*x353 - x32*x353 + x356*P[26] + x357*P[38] + x358*P[62] + x359*P[74] + x360*P[85] + x361*P[95] + x362*P[143] + x363*P[49];
   Pnew[50] = -x19*x364 - x20*x364 + x367*P[17] + x368*P[29] + x369*P[40] + x370*P[51] + x371*P[52] + x372*P[53] + x373*P[62] + x374*P[50];
   Pnew[51] = -x22*x364 - x23*x364 + x367*P[18] + x368*P[30] + x369*P[41] + x370*P[63] + x371*P[64] + x372*P[65] + x373*P[74] + x374*P[51];
   Pnew[52] = -x25*x364 - x26*x364 + x367*P[19] + x368*P[31] + x369*P[42] + x370*P[64] + x371*P[75] + x372*P[76] + x373*P[85] + x374*P[52];
   Pnew[53] = -x28*x364 - x29*x364 + x367*P[20] + x368*P[32] + x369*P[43] + x370*P[65] + x371*P[76] + x372*P[86] + x373*P[95] + x374*P[53];
   Pnew[54] = x367*P[21] + x368*P[33] + x369*P[44] + x370*P[66] + x371*P[77] + x372*P[87] + x373*P[101] + x374*P[54];
   Pnew[55] = x367*P[22] + x368*P[34] + x369*P[45] + x370*P[67] + x371*P[78] + x372*P[88] + x373*P[105] + x374*P[55];
   Pnew[56] = x367*P[23] + x368*P[35] + x369*P[46] + x370*P[68] + x371*P[79] + x372*P[89] + x373*P[110] + x374*P[56];
   Pnew[57] = -x35*x364 - x36*x364 + x367*P[24] + x368*P[36] + x369*P[47] + x370*P[69] + x371*P[80] + x372*P[90] + x373*P[125] + x374*P[57];
   Pnew[58] = x370*P[70] + x371*P[81] + x372*P[91] + x373*P[136] + x374*P[58];
   Pnew[59] = x370*P[71] + x371*P[82] + x372*P[92] + x373*P[138] + x374*P[59];
   Pnew[60] = x370*P[72] + x371*P[83] + x372*P[93] + x373*P[140] + x374*P[60];
   Pnew[61] = -x364*x38 - x364*x39 + x367*P[25] + x368*P[37] + x369*P[48] + x370*P[73] + x371*P[84] + x372*P[94] + x373*P[142] + x374*P[61];
   Pnew[62] = -x31*x364 - x32*x364 + x367*P[26] + x368*P[38] + x369*P[49] + x370*P[74] + x371*P[85] + x372*P[95] + x373*P[143] + x374*P[62];
   Pnew[63] = -x22*x375 - x23*x375 + x378*P[18] + x379*P[30] + x380*P[41] + x381*P[51] + x382*P[64] + x383*P[65] + x384*P[74] + x385*P[63];
   Pnew[64] = -x25*x375 - x26*x375 + x378*P[19] + x379*P[31] + x380*P[42] + x381*P[52] + x382*P[75] + x383*P[76] + x384*P[85] + x385*P[64];
   Pnew[65] = -x28*x375 - x29*x375 + x378*P[20] + x379*P[32] + x380*P[43] + x381*P[53] + x382*P[76] + x383*P[86] + x384*P[95] + x385*P[65];
   Pnew[66] = x378*P[21] + x379*P[33] + x380*P[44] + x381*P[54] + x382*P[77] + x383*P[87] + x384*P[101] + x385*P[66];
   Pnew[67] = x378*P[22] + x379*P[34] + x380*P[45] + x381*P[55] + x382*P[78] + x383*P[88] + x384*P[105] + x385*P[67];
   Pnew[68] = x378*P[23] + x379*P[35] + x380*P[46] + x381*P[56] + x382*P[79] + x383*P[89] + x384*P[110] + x385*P[68];
   Pnew[69] = -x35*x375 - x36*x375 + x378*P[24] + x379*P[36] + x380*P[47] + x381*P[57] + x382*P[80] + x383*P[90] + x384*P[125] + x385*P[69];
   Pnew[70] = x381*P[58] + x382*P[81] + x383*P[91] + x384*P[136] + x385*P[70];
   Pnew[71] = x381*P[59] + x382*P[82] + x383*P[92] + x384*P[138] + x385*P[71];
   Pnew[72] = x381*P[60] + x382*P[83] + x383*P[93] + x384*P[140] + x385*P[72];
   Pnew[73] = -x375*x38 - x375*x39 + x378*P[25] + x379*P[37] + x380*P[48] + x381*P[61] + x382*P[84] + x383*P[94] + x384*P[142] + x385*P[73];
   Pnew[74] = -x31*x375 - x32*x375 + x378*P[26] + x379*P[38] + x380*P[49] + x381*P[62] + x382*P[85] + x383*P[95] + x384*P[143] + x385*P[74];
   Pnew[75] = -x25*x386 - x26*x386 + x389*P[19] + x390*P[31] + x391*P[42] + x392*P[52] + x393*P[64] + x394*P[76] + x395*P[85] + x396*P[75];
   Pnew[76] = -x28*x386 - x29*x386 + x389*P[20] + x390*P[32] + x391*P[43] + x392*P[53] + x393*P[65] + x394*P[86] + x395*P[95] + x396*P[76];
   Pnew[77] = x389*P[21] + x390*P[33] + x391*P[44] + x392*P[54] + x393*P[66] + x394*P[87] + x395*P[101] + x396*P[77];
   Pnew[78] = x389*P[22] + x390*P[34] + x391*P[45] + x392*P[55] + x393*P[67] + x394*P[88] + x395*P[105] + x396*P[78];
   Pnew[79] = x389*P[23] + x390*P[35] + x391*P[46] + x392*P[56] + x393*P[68] + x394*P[89] + x395*P[110] + x396*P[79];
   Pnew[80] = -x35*x386 - x36*x386 + x389*P[24] + x390*P[36] + x391*P[47] + x392*P[57] + x393*P[69] + x394*P[90] + x395*P[125] + x396*P[80];
   Pnew[81] = x392*P[58] + x393*P[70] + x394*P[91] + x395*P[136] + x396*P[81];
   Pnew[82] = x392*P[59] + x393*P[71] + x394*P[92] + x395*P[138] + x396*P[82];
   Pnew[83] = x392*P[60] + x393*P[72] + x394*P[93] + x395*P[140] + x396*P[83];
   Pnew[84] = -x38*x386 - x386*x39 + x389*P[25] + x390*P[37] + x391*P[48] + x392*P[61] + x393*P[73] + x394*P[94] + x395*P[142] + x396*P[84];
   Pnew[85] = -x31*x386 - x32*x386 + x389*P[26] + x390*P[38] + x391*P[49] + x392*P[62] + x393*P[74] + x394*P[95] + x395*P[143] + x396*P[85];
   Pnew[86] = -x28*x397 - x29*x397 + x400*P[20] + x401*P[32] + x402*P[43] + x403*P[53] + x404*P[65] + x405*P[76] + x406*P[95] + x407*P[86];
   Pnew[87] = x400*P[21] + x401*P[33] + x402*P[44] + x403*P[54] + x404*P[66] + x405*P[77] + x406*P[101] + x407*P[87];
   Pnew[88] = x400*P[22] + x401*P[34] + x402*P[45] + x403*P[55] + x404*P[67] + x405*P[78] + x406*P[105] + x407*P[88];
   Pnew[89] = x400*P[23] + x401*P[35] + x402*P[46] + x403*P[56] + x404*P[68] + x405*P[79] + x406*P[110] + x407*P[89];
   Pnew[90] = -x35*x397 - x36*x397 + x400*P[24] + x401*P[36] + x402*P[47] + x403*P[57] + x404*P[69] + x405*P[80] + x406*P[125] + x407*P[90];
   Pnew[91] = x403*P[58] + x404*P[70] + x405*P[81] + x406*P[136] + x407*P[91];
   Pnew[92] = x403*P[59] + x404*P[71] + x405*P[82] + x406*P[138] + x407*P[92];
   Pnew[93] = x403*P[60] + x404*P[72] + x405*P[83] + x406*P[140] + x407*P[93];
   Pnew[94] = -x38*x397 - x39*x397 + x400*P[25] + x401*P[37] + x402*P[48] + x403*P[61] + x404*P[73] + x405*P[84] + x406*P[142] + x407*P[94];
   Pnew[95] = -x31*x397 - x32*x397 + x400*P[26] + x401*P[38] + x402*P[49] + x403*P[62] + x404*P[74] + x405*P[85] + x406*P[143] + x407*P[95];
   Pnew[96] = x410*P[21] + x411*P[33] + x412*P[44] + x413*P[54] + x414*P[66] + x415*P[77] + x416*P[87] + x417*P[101] + P[96];
   Pnew[97] = x417*P[114] + P[97];
   Pnew[98] = x417*P[118] + P[98];
   Pnew[99] = x417*P[128] + P[99];
   Pnew[100] = x413*P[58] + x414*P[70] + x415*P[81] + x416*P[91] + x417*P[136] + P[100];
   Pnew[101] = -x31*x418 - x32*x418 + x410*P[26] + x411*P[38] + x412*P[49] + x413*P[62] + x414*P[74] + x415*P[85] + x416*P[95] + x417*P[143] + P[101];
   Pnew[102] = x421*P[22] + x422*P[34] + x423*P[45] + x424*P[55] + x425*P[67] + x426*P[78] + x427*P[88] + x428*P[105] + P[102];
   Pnew[103] = x428*P[131] + P[103];
   Pnew[104] = x424*P[59] + x425*P[71] + x426*P[82] + x427*P[92] + x428*P[138] + P[104];
   Pnew[105] = -x31*x429 - x32*x429 + x421*P[26] + x422*P[38] + x423*P[49] + x424*P[62] + x425*P[74] + x426*P[85] + x427*P[95] + x428*P[143] + P[105];
   Pnew[106] = x432*P[23] + x433*P[35] + x434*P[46] + x435*P[56] + x436*P[68] + x437*P[79] + x438*P[89] + x439*P[110] + P[106];
   Pnew[107] = x439*P[122] + P[107];
   Pnew[108] = x439*P[134] + P[108];
   Pnew[109] = x435*P[60] + x436*P[72] + x437*P[83] + x438*P[93] + x439*P[140] + P[109];
   Pnew[110] = -x31*x440 - x32*x440 + x432*P[26] + x433*P[38] + x434*P[49] + x435*P[62] + x436*P[74] + x437*P[85] + x438*P[95] + x439*P[143] + P[110];
   Pnew[111] = x443*P[114] + P[111];
   Pnew[112] = x443*P[128] + P[112];
   Pnew[113] = x443*P[136] + x444*P[58] + x445*P[70] + x446*P[81] + x447*P[91] + P[113];
   Pnew[114] = -x31*x448 - x32*x448 + x443*P[143] + x444*P[62] + x445*P[74] + x446*P[85] + x447*P[95] + (-x441*H[9] - x442*H[1])*P[26] + (-x441*H[10] - x442*H[2])*P[38] + (-x441*H[11] - x442*H[3])*P[49] + P[114];
   Pnew[115] = x451*P[118] + P[115];
   Pnew[116] = x451*P[131] + P[116];
   Pnew[117] = x451*P[138] + x452*P[59] + x453*P[71] + x454*P[82] + x455*P[92] + P[117];
   Pnew[118] = -x31*x456 - x32*x456 + x451*P[143] + x452*P[62] + x453*P[74] + x454*P[85] + x455*P[95] + (-x449*H[9] - x450*H[1])*P[26] + (-x449*H[10] - x450*H[2])*P[38] + (-x449*H[11] - x450*H[3])*P[49] + P[118];
   Pnew[119] = x459*P[122] + P[119];
   Pnew[120] = x459*P[134] + P[120];
   Pnew[121] = x459*P[140] + x460*P[60] + x461*P[72] + x462*P[83] + x463*P[93] + P[121];
   Pnew[122] = -x31*x464 - x32*x464 + x459*P[143] + x460*P[62] + x461*P[74] + x462*P[85] + x463*P[95] + (-x457*H[9] - x458*H[1])*P[26] + (-x457*H[10] - x458*H[2])*P[38] + (-x457*H[11] - x458*H[3])*P[49] + P[122];
   Pnew[123] = -x36*x465 + x466*P[123] + x469*P[24] + x470*P[36] + x471*P[47] + x472*P[57] + x473*P[69] + x474*P[80] + x475*P[90] + x476*P[125];
   Pnew[124] = -x39*x465 + x466*P[124] + x469*P[25] + x470*P[37] + x471*P[48] + x472*P[61] + x473*P[73] + x474*P[84] + x475*P[94] + x476*P[142];
   Pnew[125] = -x32*x465 + x466*P[125] + x469*P[26] + x470*P[38] + x471*P[49] + x472*P[62] + x473*P[74] + x474*P[85] + x475*P[95] + x476*P[143];
   Pnew[126] = x479*P[128] + P[126];
   Pnew[127] = x479*P[136] + x480*P[58] + x481*P[70] + x482*P[81] + x483*P[91] + P[127];
   Pnew[128] = -x31*x484 - x32*x484 + x479*P[143] + x480*P[62] + x481*P[74] + x482*P[85] + x483*P[95] + (-x477*H[9] - x478*H[1])*P[26] + (-x477*H[10] - x478*H[2])*P[38] + (-x477*H[11] - x478*H[3])*P[49] + P[128];
   Pnew[129] = x487*P[131] + P[129];
   Pnew[130] = x487*P[138] + x488*P[59] + x489*P[71] + x490*P[82] + x491*P[92] + P[130];
   Pnew[131] = -x31*x492 - x32*x492 + x487*P[143] + x488*P[62] + x489*P[74] + x490*P[85] + x491*P[95] + (-x485*H[9] - x486*H[1])*P[26] + (-x485*H[10] - x486*H[2])*P[38] + (-x485*H[11] - x486*H[3])*P[49] + P[131];
   Pnew[132] = x495*P[134] + P[132];
   Pnew[133] = x495*P[140] + x496*P[60] + x497*P[72] + x498*P[83] + x499*P[93] + P[133];
   Pnew[134] = -x31*x500 - x32*x500 + x495*P[143] + x496*P[62] + x497*P[74] + x498*P[85] + x499*P[95] + (-x493*H[9] - x494*H[1])*P[26] + (-x493*H[10] - x494*H[2])*P[38] + (-x493*H[11] - x494*H[3])*P[49] + P[134];
   Pnew[135] = x503*P[58] + x504*P[70] + x505*P[81] + x506*P[91] + x507*P[136] + P[135];
   Pnew[136] = -x31*x508 - x32*x508 + x503*P[62] + x504*P[74] + x505*P[85] + x506*P[95] + x507*P[143] + (-x501*H[9] - x502*H[1])*P[26] + (-x501*H[10] - x502*H[2])*P[38] + (-x501*H[11] - x502*H[3])*P[49] + P[136];
   Pnew[137] = x511*P[59] + x512*P[71] + x513*P[82] + x514*P[92] + x515*P[138] + P[137];
   Pnew[138] = -x31*x516 - x32*x516 + x511*P[62] + x512*P[74] + x513*P[85] + x514*P[95] + x515*P[143] + (-x509*H[9] - x510*H[1])*P[26] + (-x509*H[10] - x510*H[2])*P[38] + (-x509*H[11] - x510*H[3])*P[49] + P[138];
   Pnew[139] = x519*P[60] + x520*P[72] + x521*P[83] + x522*P[93] + x523*P[140] + P[139];
   Pnew[140] = -x31*x524 - x32*x524 + x519*P[62] + x520*P[74] + x521*P[85] + x522*P[95] + x523*P[143] + (-x517*H[9] - x518*H[1])*P[26] + (-x517*H[10] - x518*H[2])*P[38] + (-x517*H[11] - x518*H[3])*P[49] + P[140];
   Pnew[141] = -x38*x525 + x526*P[141] + x529*P[25] + x530*P[37] + x531*P[48] + x532*P[61] + x533*P[73] + x534*P[84] + x535*P[94] + x536*P[142];
   Pnew[142] = -x31*x525 + x526*P[142] + x529*P[26] + x530*P[38] + x531*P[49] + x532*P[62] + x533*P[74] + x534*P[85] + x535*P[95] + x536*P[143];
   Pnew[143] = -x31*x537 - x32*x537 + (-x538*H[9] - x539*H[1])*P[26] + (-x538*H[10] - x539*H[2])*P[38] + (-x538*H[11] - x539*H[3])*P[49] + (-x538*H[12] - x539*H[4])*P[62] + (-x538*H[13] - x539*H[5])*P[74] + (-x538*H[14] - x539*H[6])*P[85] + (-x538*H[15] - x539*H[7])*P[95] + (-x538*H[16] - x539*H[8] + 1)*P[143];
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
