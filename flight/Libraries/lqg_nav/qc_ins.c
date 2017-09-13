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
#define NUMX 18  // number of state variables
#define NUMF 60  // number of non-zero variables stored in the F matrix
#define NUMH 31  // number of non-zero variables stored in the H matrix
#define NUMP 87  // number of non-zero variables in the upper triangular of covariance
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
		float beta_t;
		float tau;
		float mu;
		float Be[3];
	} __attribute__((__packed__)) params;

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
	Q[15] = Q[16] = Q[17] = 1e-6f;      // Bias states

	// Store the observation noises
	qcins_state->R[0] = 100;
	qcins_state->R[1] = qcins_state->R[2] = qcins_state->R[3] = 1e6;
	qcins_state->R[4] = qcins_state->R[5] = qcins_state->R[6] = 1e3;
	qcins_state->R[7] = qcins_state->R[8] = qcins_state->R[9] = 1e3;

	// Defaults for the parameters
	qcins_state->params.g = 9.81f;
	qcins_state->params.beta_r = 10000.0f * DEG2RAD;
	qcins_state->params.beta_p = 10000.0f * DEG2RAD;
	qcins_state->params.beta_y1 = 1000.0f * DEG2RAD;
	qcins_state->params.beta_y2 = 1000.0f * DEG2RAD;
	qcins_state->params.beta_t = 2.0f;
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

bool qcins_set_thrust(uintptr_t qcins_handle, const float beta_t_new)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->params.beta_t = beta_t_new;
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

bool qcins_predict(uintptr_t qcins_handle, const float roll, const float pitch, const float yaw, const float throttle, float Ts)
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	const float u[4] = {roll, pitch, yaw, throttle};

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

	rate[0] = qcins_state->x[8];
	rate[1] = qcins_state->x[9];
	rate[2] = qcins_state->x[10];
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

void state_prediction(const float * restrict x, const float *u, float Ts, const float *restrict param, float * restrict xnew)
{
   const float x0 = 2*param[0]*param[5]*x[14];
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
   const float x17 = x16*param[7];
   const float x18 = -x12;
   const float x19 = x10 + x14 + x18 + x9;
   const float x20 = x3 + 2*x[5]*x[6];
   const float x21 = x1 - x2;
   const float x22 = x19*x[1] + x20*x[2] - x21*x5;
   const float x23 = x22*param[7];
   const float x24 = 2*param[7];
   const float x25 = (1.0F/2.0F)*x[5];
   const float x26 = (1.0F/2.0F)*x[6];
   const float x27 = (1.0F/2.0F)*x[7];
   const float x28 = (1.0F/2.0F)*x[4];
   const float x29 = Ts/param[6];
   xnew[0] = Ts*x[3] + x[0];
   xnew[1] = Ts*(-x0*(x1 + x2) + x17*x4 - x19*x23) + x[1];
   xnew[2] = Ts*(x0*(x6 - x7) - x15*x17 - x20*x23) + x[2];
   xnew[3] = Ts*(-x16*x24*x8 + x21*x22*x24 - (x11 + x13 + x18)*param[0]*param[5]*x[14] + param[0]) + x[3];
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
   const float x12 = x4*param[7];
   const float x13 = x5 - x6;
   const float x14 = x10 + x13 + x7;
   const float x15 = x1 + x3;
   const float x16 = x15*param[7];
   const float x17 = Ts*(-x11*x16 - x12*x14);
   const float x18 = x0*x[5];
   const float x19 = 2*x[6];
   const float x20 = x19*x[7];
   const float x21 = x18 + x20;
   const float x22 = x21*param[7];
   const float x23 = x0*x[6];
   const float x24 = x2*x[7];
   const float x25 = -x23 + x24;
   const float x26 = x25*param[7];
   const float x27 = Ts*(-x11*x26 - x22*x4);
   const float x28 = 2*x[1];
   const float x29 = 2*x[2];
   const float x30 = 2*x[3];
   const float x31 = -x28*x[7] + x29*x[4] + x30*x[5];
   const float x32 = x28*x[4];
   const float x33 = x29*x[7];
   const float x34 = x30*x[6];
   const float x35 = x32 + x33 - x34;
   const float x36 = x11*param[7];
   const float x37 = param[0]*param[5]*x[14];
   const float x38 = x19*x37;
   const float x39 = 2*x[7];
   const float x40 = (x14*x[2] + x21*x[3] + x4*x[1])*param[7];
   const float x41 = x39*x40;
   const float x42 = (x11*x[1] + x15*x[2] + x25*x[3])*param[7];
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
   const float x60 = Ts*param[0]*param[5];
   const float x61 = Ts*(-x14*x22 - x15*x26);
   const float x62 = x14*param[7];
   const float x63 = x56 - x58 - x59;
   const float x64 = (1.0F/2.0F)*Ts;
   const float x65 = x64*x[8];
   const float x66 = -x65;
   const float x67 = x64*x[9];
   const float x68 = -x67;
   const float x69 = x64*x[10];
   const float x70 = -x69;
   const float x71 = x64*x[5];
   const float x72 = -x71;
   const float x73 = x64*x[6];
   const float x74 = -x73;
   const float x75 = x64*x[7];
   const float x76 = -x75;
   const float x77 = x64*x[4];
   const float x78 = -Ts/param[6];
   const float x79 = -x49*param[7];
   const float x80 = -x31*param[7];
   const float x81 = 2*param[8];
   const float x82 = x81*x[4];
   const float x83 = 2*param[9];
   const float x84 = x83*x[7];
   const float x85 = 2*param[10];
   const float x86 = x85*x[6];
   const float x87 = x82 + x84 - x86;
   const float x88 = x81*x[5] + x83*x[6] + x85*x[7];
   const float x89 = x81*x[6];
   const float x90 = x83*x[5];
   const float x91 = x85*x[4];
   const float x92 = x81*x[7];
   const float x93 = x83*x[4];
   const float x94 = x85*x[5];
   const float x95 = -x92 + x93 + x94;
   const float x96 = x89 - x90 + x91;
   F[0] = Ts;
   F[1] = Ts*(-(x11*x11)*param[7] - (x4*x4)*param[7]);
   F[2] = x17;
   F[3] = x27;
   F[4] = Ts*(-x12*x31 - x35*x36 + x44);
   F[5] = Ts*(-x12*x48 - x36*x49 + x50);
   F[6] = Ts*(-x12*x49 - x36*x51 + x55);
   F[7] = Ts*(-x12*x57 - x31*x36 - x56 + x58 + x59);
   F[8] = -x60*(x23 + x24);
   F[9] = x17;
   F[10] = Ts*(-(x14*x14)*param[7] - (x15*x15)*param[7]);
   F[11] = x61;
   F[12] = Ts*(-x16*x35 - x31*x62 + x63);
   F[13] = Ts*(-x16*x49 - x48*x62 + x52 + x53 - x54);
   F[14] = Ts*(-x16*x51 - x49*x62 + x50);
   F[15] = Ts*(-x16*x31 + x44 - x57*x62);
   F[16] = -x60*(-x18 + x20);
   F[17] = x27;
   F[18] = x61;
   F[19] = Ts*(-(x21*x21)*param[7] - (x25*x25)*param[7]);
   F[20] = Ts*(-x22*x31 - x26*x35 + x55);
   F[21] = Ts*(-x22*x48 - x26*x49 + x63);
   F[22] = Ts*(-x22*x49 - x26*x51 + x38 - x41 + x43);
   F[23] = Ts*(-x22*x57 - x26*x31 + x50);
   F[24] = -x60*(x13 + x8 + x9);
   F[25] = x66;
   F[26] = x68;
   F[27] = x70;
   F[28] = x72;
   F[29] = x74;
   F[30] = x76;
   F[31] = x65;
   F[32] = x69;
   F[33] = x68;
   F[34] = x77;
   F[35] = x76;
   F[36] = x73;
   F[37] = x67;
   F[38] = x70;
   F[39] = x65;
   F[40] = x75;
   F[41] = x77;
   F[42] = x72;
   F[43] = x69;
   F[44] = x67;
   F[45] = x66;
   F[46] = x74;
   F[47] = x71;
   F[48] = x77;
   F[49] = Ts*param[1];
   F[50] = Ts*param[2];
   F[51] = Ts*(param[3] - param[4]);
   F[52] = -Ts*param[4];
   F[53] = x78;
   F[54] = x78;
   F[55] = x78;
   F[56] = x78;
   F[57] = x78;
   F[58] = x78;
   F[59] = x78;
   H[0] = 1;
   H[1] = -x36;
   H[2] = -x16;
   H[3] = -x26;
   H[4] = -x35*param[7];
   H[5] = x79;
   H[6] = -x51*param[7];
   H[7] = x80;
   H[8] = -x12;
   H[9] = -x62;
   H[10] = -x22;
   H[11] = x80;
   H[12] = -x48*param[7];
   H[13] = x79;
   H[14] = -x57*param[7];
   H[15] = -param[0]*param[5];
   H[16] = 1;
   H[17] = 1;
   H[18] = 1;
   H[19] = x87;
   H[20] = x88;
   H[21] = -x89 + x90 - x91;
   H[22] = x95;
   H[23] = x95;
   H[24] = x96;
   H[25] = x88;
   H[26] = -x82 - x84 + x86;
   H[27] = x96;
   H[28] = x92 - x93 - x94;
   H[29] = x87;
   H[30] = x88;
}

// TODO: rework F definition so that Ts is already multiplied into all the terms
// to reduce the memory requirement (currently duplicates F)
void covariance_prediction(const float *restrict P, const float * restrict F, const float * restrict Q, float * restrict Pnew)
{
   const float x0 = F[0]*P[33] + P[3];
   const float x1 = F[0]*P[24] + P[2];
   const float x2 = F[0]*P[34] + P[4];
   const float x3 = F[0]*P[35] + P[5];
   const float x4 = F[0]*P[36] + P[6];
   const float x5 = F[0]*P[37] + P[7];
   const float x6 = F[0]*P[41] + P[11];
   const float x7 = F[1] + 1;
   const float x8 = F[0]*P[14] + P[1];
   const float x9 = F[10] + 1;
   const float x10 = F[19] + 1;
   const float x11 = F[0]*P[38] + P[8];
   const float x12 = F[0]*P[39] + P[9];
   const float x13 = F[0]*P[40] + P[10];
   const float x14 = F[59] + 1;
   const float x15 = x7*P[13] + F[2]*P[23] + F[3]*P[24] + F[4]*P[25] + F[5]*P[26] + F[6]*P[27] + F[7]*P[28] + F[8]*P[32];
   const float x16 = x7*P[14] + F[2]*P[24] + F[3]*P[33] + F[4]*P[34] + F[5]*P[35] + F[6]*P[36] + F[7]*P[37] + F[8]*P[41];
   const float x17 = x7*P[15] + F[2]*P[25] + F[3]*P[34] + F[4]*P[42] + F[5]*P[43] + F[6]*P[44] + F[7]*P[45] + F[8]*P[49];
   const float x18 = x7*P[16] + F[2]*P[26] + F[3]*P[35] + F[4]*P[43] + F[5]*P[50] + F[6]*P[51] + F[7]*P[52] + F[8]*P[56];
   const float x19 = x7*P[17] + F[2]*P[27] + F[3]*P[36] + F[4]*P[44] + F[5]*P[51] + F[6]*P[57] + F[7]*P[58] + F[8]*P[62];
   const float x20 = x7*P[18] + F[2]*P[28] + F[3]*P[37] + F[4]*P[45] + F[5]*P[52] + F[6]*P[58] + F[7]*P[63] + F[8]*P[67];
   const float x21 = x7*P[22] + F[2]*P[32] + F[3]*P[41] + F[4]*P[49] + F[5]*P[56] + F[6]*P[62] + F[7]*P[67] + F[8]*P[83];
   const float x22 = x7*P[12] + F[2]*P[13] + F[3]*P[14] + F[4]*P[15] + F[5]*P[16] + F[6]*P[17] + F[7]*P[18] + F[8]*P[22];
   const float x23 = x7*P[19] + F[2]*P[29] + F[3]*P[38] + F[4]*P[46] + F[5]*P[53] + F[6]*P[59] + F[7]*P[64];
   const float x24 = x7*P[20] + F[2]*P[30] + F[3]*P[39] + F[4]*P[47] + F[5]*P[54] + F[6]*P[60] + F[7]*P[65];
   const float x25 = x7*P[21] + F[2]*P[31] + F[3]*P[40] + F[4]*P[48] + F[5]*P[55] + F[6]*P[61] + F[7]*P[66];
   const float x26 = x9*P[13] + F[9]*P[12] + F[11]*P[14] + F[12]*P[15] + F[13]*P[16] + F[14]*P[17] + F[15]*P[18] + F[16]*P[22];
   const float x27 = x9*P[24] + F[9]*P[14] + F[11]*P[33] + F[12]*P[34] + F[13]*P[35] + F[14]*P[36] + F[15]*P[37] + F[16]*P[41];
   const float x28 = x9*P[25] + F[9]*P[15] + F[11]*P[34] + F[12]*P[42] + F[13]*P[43] + F[14]*P[44] + F[15]*P[45] + F[16]*P[49];
   const float x29 = x9*P[26] + F[9]*P[16] + F[11]*P[35] + F[12]*P[43] + F[13]*P[50] + F[14]*P[51] + F[15]*P[52] + F[16]*P[56];
   const float x30 = x9*P[27] + F[9]*P[17] + F[11]*P[36] + F[12]*P[44] + F[13]*P[51] + F[14]*P[57] + F[15]*P[58] + F[16]*P[62];
   const float x31 = x9*P[28] + F[9]*P[18] + F[11]*P[37] + F[12]*P[45] + F[13]*P[52] + F[14]*P[58] + F[15]*P[63] + F[16]*P[67];
   const float x32 = x9*P[32] + F[9]*P[22] + F[11]*P[41] + F[12]*P[49] + F[13]*P[56] + F[14]*P[62] + F[15]*P[67] + F[16]*P[83];
   const float x33 = x9*P[23] + F[9]*P[13] + F[11]*P[24] + F[12]*P[25] + F[13]*P[26] + F[14]*P[27] + F[15]*P[28] + F[16]*P[32];
   const float x34 = x9*P[29] + F[9]*P[19] + F[11]*P[38] + F[12]*P[46] + F[13]*P[53] + F[14]*P[59] + F[15]*P[64];
   const float x35 = x9*P[30] + F[9]*P[20] + F[11]*P[39] + F[12]*P[47] + F[13]*P[54] + F[14]*P[60] + F[15]*P[65];
   const float x36 = x9*P[31] + F[9]*P[21] + F[11]*P[40] + F[12]*P[48] + F[13]*P[55] + F[14]*P[61] + F[15]*P[66];
   const float x37 = x10*P[34] + F[17]*P[15] + F[18]*P[25] + F[20]*P[42] + F[21]*P[43] + F[22]*P[44] + F[23]*P[45] + F[24]*P[49];
   const float x38 = x10*P[35] + F[17]*P[16] + F[18]*P[26] + F[20]*P[43] + F[21]*P[50] + F[22]*P[51] + F[23]*P[52] + F[24]*P[56];
   const float x39 = x10*P[36] + F[17]*P[17] + F[18]*P[27] + F[20]*P[44] + F[21]*P[51] + F[22]*P[57] + F[23]*P[58] + F[24]*P[62];
   const float x40 = x10*P[37] + F[17]*P[18] + F[18]*P[28] + F[20]*P[45] + F[21]*P[52] + F[22]*P[58] + F[23]*P[63] + F[24]*P[67];
   const float x41 = x10*P[41] + F[17]*P[22] + F[18]*P[32] + F[20]*P[49] + F[21]*P[56] + F[22]*P[62] + F[23]*P[67] + F[24]*P[83];
   const float x42 = x10*P[38] + F[17]*P[19] + F[18]*P[29] + F[20]*P[46] + F[21]*P[53] + F[22]*P[59] + F[23]*P[64];
   const float x43 = x10*P[39] + F[17]*P[20] + F[18]*P[30] + F[20]*P[47] + F[21]*P[54] + F[22]*P[60] + F[23]*P[65];
   const float x44 = x10*P[40] + F[17]*P[21] + F[18]*P[31] + F[20]*P[48] + F[21]*P[55] + F[22]*P[61] + F[23]*P[66];
   const float x45 = F[25]*P[43] + F[26]*P[44] + F[27]*P[45] + F[28]*P[46] + F[29]*P[47] + F[30]*P[48] + P[42];
   const float x46 = F[25]*P[53] + F[26]*P[59] + F[27]*P[64] + F[28]*P[68] + P[46];
   const float x47 = F[25]*P[54] + F[26]*P[60] + F[27]*P[65] + F[29]*P[72] + P[47];
   const float x48 = F[25]*P[55] + F[26]*P[61] + F[27]*P[66] + F[30]*P[74] + P[48];
   const float x49 = F[25]*P[50] + F[26]*P[51] + F[27]*P[52] + F[28]*P[53] + F[29]*P[54] + F[30]*P[55] + P[43];
   const float x50 = F[25]*P[51] + F[26]*P[57] + F[27]*P[58] + F[28]*P[59] + F[29]*P[60] + F[30]*P[61] + P[44];
   const float x51 = F[25]*P[52] + F[26]*P[58] + F[27]*P[63] + F[28]*P[64] + F[29]*P[65] + F[30]*P[66] + P[45];
   const float x52 = F[49]*P[69];
   const float x53 = F[50]*P[70];
   const float x54 = F[51]*P[75];
   const float x55 = F[52]*P[76];
   const float x56 = F[31]*P[43] + F[32]*P[51] + F[33]*P[52] + F[34]*P[53] + F[35]*P[54] + F[36]*P[55] + P[50];
   const float x57 = F[31]*P[46] + F[32]*P[59] + F[33]*P[64] + F[34]*P[68] + P[53];
   const float x58 = F[31]*P[47] + F[32]*P[60] + F[33]*P[65] + F[35]*P[72] + P[54];
   const float x59 = F[31]*P[48] + F[32]*P[61] + F[33]*P[66] + F[36]*P[74] + P[55];
   const float x60 = F[31]*P[42] + F[32]*P[44] + F[33]*P[45] + F[34]*P[46] + F[35]*P[47] + F[36]*P[48] + P[43];
   const float x61 = F[31]*P[44] + F[32]*P[57] + F[33]*P[58] + F[34]*P[59] + F[35]*P[60] + F[36]*P[61] + P[51];
   const float x62 = F[31]*P[45] + F[32]*P[58] + F[33]*P[63] + F[34]*P[64] + F[35]*P[65] + F[36]*P[66] + P[52];
   const float x63 = F[37]*P[44] + F[38]*P[51] + F[39]*P[58] + F[40]*P[59] + F[41]*P[60] + F[42]*P[61] + P[57];
   const float x64 = F[37]*P[46] + F[38]*P[53] + F[39]*P[64] + F[40]*P[68] + P[59];
   const float x65 = F[37]*P[47] + F[38]*P[54] + F[39]*P[65] + F[41]*P[72] + P[60];
   const float x66 = F[37]*P[48] + F[38]*P[55] + F[39]*P[66] + F[42]*P[74] + P[61];
   const float x67 = F[37]*P[42] + F[38]*P[43] + F[39]*P[45] + F[40]*P[46] + F[41]*P[47] + F[42]*P[48] + P[44];
   const float x68 = F[37]*P[43] + F[38]*P[50] + F[39]*P[52] + F[40]*P[53] + F[41]*P[54] + F[42]*P[55] + P[51];
   const float x69 = F[37]*P[45] + F[38]*P[52] + F[39]*P[63] + F[40]*P[64] + F[41]*P[65] + F[42]*P[66] + P[58];
   const float x70 = F[43]*P[46] + F[44]*P[53] + F[45]*P[59] + F[46]*P[68] + P[64];
   const float x71 = F[43]*P[47] + F[44]*P[54] + F[45]*P[60] + F[47]*P[72] + P[65];
   const float x72 = F[43]*P[48] + F[44]*P[55] + F[45]*P[61] + F[48]*P[74] + P[66];
   const float x73 = F[49]*P[77] + P[69];
   const float x74 = F[49]*P[78] + P[71];
   const float x75 = F[53] + 1;
   const float x76 = F[55] + 1;
   const float x77 = F[51]*P[81] + F[52]*P[82] + P[75];
   const float x78 = F[51]*P[82] + F[52]*P[86] + P[76];
   const float x79 = F[57] + 1;
   const float x80 = x75*P[78] + F[54]*P[84];
   const float x81 = x76*P[80] + F[56]*P[85];
   const float x82 = x79*P[82] + F[58]*P[86];
   Pnew[0] = x0*F[0] + F[0]*P[3] + P[0] + Q[0];
   Pnew[1] = x0*F[3] + x1*F[2] + x2*F[4] + x3*F[5] + x4*F[6] + x5*F[7] + x6*F[8] + x7*x8;
   Pnew[2] = x0*F[11] + x1*x9 + x2*F[12] + x3*F[13] + x4*F[14] + x5*F[15] + x6*F[16] + x8*F[9];
   Pnew[3] = x0*x10 + x1*F[18] + x2*F[20] + x3*F[21] + x4*F[22] + x5*F[23] + x6*F[24] + x8*F[17];
   Pnew[4] = x11*F[28] + x12*F[29] + x13*F[30] + x2 + x3*F[25] + x4*F[26] + x5*F[27];
   Pnew[5] = x11*F[34] + x12*F[35] + x13*F[36] + x2*F[31] + x3 + x4*F[32] + x5*F[33];
   Pnew[6] = x11*F[40] + x12*F[41] + x13*F[42] + x2*F[37] + x3*F[38] + x4 + x5*F[39];
   Pnew[7] = x11*F[46] + x12*F[47] + x13*F[48] + x2*F[43] + x3*F[44] + x4*F[45] + x5;
   Pnew[8] = x11;
   Pnew[9] = x12;
   Pnew[10] = x13;
   Pnew[11] = x14*x6;
   Pnew[12] = x15*F[2] + x16*F[3] + x17*F[4] + x18*F[5] + x19*F[6] + x20*F[7] + x21*F[8] + x22*x7 + Q[1];
   Pnew[13] = x15*x9 + x16*F[11] + x17*F[12] + x18*F[13] + x19*F[14] + x20*F[15] + x21*F[16] + x22*F[9];
   Pnew[14] = x10*x16 + x15*F[18] + x17*F[20] + x18*F[21] + x19*F[22] + x20*F[23] + x21*F[24] + x22*F[17];
   Pnew[15] = x17 + x18*F[25] + x19*F[26] + x20*F[27] + x23*F[28] + x24*F[29] + x25*F[30];
   Pnew[16] = x17*F[31] + x18 + x19*F[32] + x20*F[33] + x23*F[34] + x24*F[35] + x25*F[36];
   Pnew[17] = x17*F[37] + x18*F[38] + x19 + x20*F[39] + x23*F[40] + x24*F[41] + x25*F[42];
   Pnew[18] = x17*F[43] + x18*F[44] + x19*F[45] + x20 + x23*F[46] + x24*F[47] + x25*F[48];
   Pnew[19] = x23;
   Pnew[20] = x24;
   Pnew[21] = x25;
   Pnew[22] = x14*x21;
   Pnew[23] = x26*F[9] + x27*F[11] + x28*F[12] + x29*F[13] + x30*F[14] + x31*F[15] + x32*F[16] + x33*x9 + Q[2];
   Pnew[24] = x10*x27 + x26*F[17] + x28*F[20] + x29*F[21] + x30*F[22] + x31*F[23] + x32*F[24] + x33*F[18];
   Pnew[25] = x28 + x29*F[25] + x30*F[26] + x31*F[27] + x34*F[28] + x35*F[29] + x36*F[30];
   Pnew[26] = x28*F[31] + x29 + x30*F[32] + x31*F[33] + x34*F[34] + x35*F[35] + x36*F[36];
   Pnew[27] = x28*F[37] + x29*F[38] + x30 + x31*F[39] + x34*F[40] + x35*F[41] + x36*F[42];
   Pnew[28] = x28*F[43] + x29*F[44] + x30*F[45] + x31 + x34*F[46] + x35*F[47] + x36*F[48];
   Pnew[29] = x34;
   Pnew[30] = x35;
   Pnew[31] = x36;
   Pnew[32] = x14*x32;
   Pnew[33] = x10*(x10*P[33] + F[17]*P[14] + F[18]*P[24] + F[20]*P[34] + F[21]*P[35] + F[22]*P[36] + F[23]*P[37] + F[24]*P[41]) + x37*F[20] + x38*F[21] + x39*F[22] + x40*F[23] + x41*F[24] + (x10*P[14] + F[17]*P[12] + F[18]*P[13] + F[20]*P[15] + F[21]*P[16] + F[22]*P[17] + F[23]*P[18] + F[24]*P[22])*F[17] + (x10*P[24] + F[17]*P[13] + F[18]*P[23] + F[20]*P[25] + F[21]*P[26] + F[22]*P[27] + F[23]*P[28] + F[24]*P[32])*F[18] + Q[3];
   Pnew[34] = x37 + x38*F[25] + x39*F[26] + x40*F[27] + x42*F[28] + x43*F[29] + x44*F[30];
   Pnew[35] = x37*F[31] + x38 + x39*F[32] + x40*F[33] + x42*F[34] + x43*F[35] + x44*F[36];
   Pnew[36] = x37*F[37] + x38*F[38] + x39 + x40*F[39] + x42*F[40] + x43*F[41] + x44*F[42];
   Pnew[37] = x37*F[43] + x38*F[44] + x39*F[45] + x40 + x42*F[46] + x43*F[47] + x44*F[48];
   Pnew[38] = x42;
   Pnew[39] = x43;
   Pnew[40] = x44;
   Pnew[41] = x14*x41;
   Pnew[42] = x45 + x46*F[28] + x47*F[29] + x48*F[30] + x49*F[25] + x50*F[26] + x51*F[27] + Q[4];
   Pnew[43] = x45*F[31] + x46*F[34] + x47*F[35] + x48*F[36] + x49 + x50*F[32] + x51*F[33];
   Pnew[44] = x45*F[37] + x46*F[40] + x47*F[41] + x48*F[42] + x49*F[38] + x50 + x51*F[39];
   Pnew[45] = x45*F[43] + x46*F[46] + x47*F[47] + x48*F[48] + x49*F[44] + x50*F[45] + x51;
   Pnew[46] = x46 + x52*F[28];
   Pnew[47] = x47 + x53*F[28];
   Pnew[48] = x48 + x54*F[30] + x55*F[30];
   Pnew[49] = x14*(F[25]*P[56] + F[26]*P[62] + F[27]*P[67] + P[49]);
   Pnew[50] = x56 + x57*F[34] + x58*F[35] + x59*F[36] + x60*F[31] + x61*F[32] + x62*F[33] + Q[5];
   Pnew[51] = x56*F[38] + x57*F[40] + x58*F[41] + x59*F[42] + x60*F[37] + x61 + x62*F[39];
   Pnew[52] = x56*F[44] + x57*F[46] + x58*F[47] + x59*F[48] + x60*F[43] + x61*F[45] + x62;
   Pnew[53] = x52*F[34] + x57;
   Pnew[54] = x53*F[34] + x58;
   Pnew[55] = x54*F[36] + x55*F[36] + x59;
   Pnew[56] = x14*(F[31]*P[49] + F[32]*P[62] + F[33]*P[67] + P[56]);
   Pnew[57] = x63 + x64*F[40] + x65*F[41] + x66*F[42] + x67*F[37] + x68*F[38] + x69*F[39] + Q[6];
   Pnew[58] = x63*F[45] + x64*F[46] + x65*F[47] + x66*F[48] + x67*F[43] + x68*F[44] + x69;
   Pnew[59] = x52*F[40] + x64;
   Pnew[60] = x53*F[40] + x65;
   Pnew[61] = x54*F[42] + x55*F[42] + x66;
   Pnew[62] = x14*(F[37]*P[49] + F[38]*P[56] + F[39]*P[67] + P[62]);
   Pnew[63] = x70*F[46] + x71*F[47] + x72*F[48] + (F[43]*P[42] + F[44]*P[43] + F[45]*P[44] + F[46]*P[46] + F[47]*P[47] + F[48]*P[48] + P[45])*F[43] + (F[43]*P[43] + F[44]*P[50] + F[45]*P[51] + F[46]*P[53] + F[47]*P[54] + F[48]*P[55] + P[52])*F[44] + (F[43]*P[44] + F[44]*P[51] + F[45]*P[57] + F[46]*P[59] + F[47]*P[60] + F[48]*P[61] + P[58])*F[45] + F[43]*P[45] + F[44]*P[52] + F[45]*P[58] + F[46]*P[64] + F[47]*P[65] + F[48]*P[66] + P[63] + Q[7];
   Pnew[64] = x52*F[46] + x70;
   Pnew[65] = x53*F[46] + x71;
   Pnew[66] = x54*F[48] + x55*F[48] + x72;
   Pnew[67] = x14*(F[43]*P[49] + F[44]*P[56] + F[45]*P[62] + P[67]);
   Pnew[68] = x52 + x73*F[49] + P[68] + Q[8];
   Pnew[69] = x73*x75 + x74*F[54];
   Pnew[70] = x76*P[70];
   Pnew[71] = x74;
   Pnew[72] = (F[50]*F[50])*P[79] + P[72] + Q[9];
   Pnew[73] = F[50]*P[80] + P[73];
   Pnew[74] = x54 + x55 + x77*F[51] + x78*F[52] + P[74] + Q[10];
   Pnew[75] = x77*x79 + x78*F[58];
   Pnew[76] = x78;
   Pnew[77] = x75*(x75*P[77] + F[54]*P[78]) + x80*F[54] + Q[11];
   Pnew[78] = x80;
   Pnew[79] = x76*(x76*P[79] + F[56]*P[80]) + x81*F[56] + Q[12];
   Pnew[80] = x81;
   Pnew[81] = x79*(x79*P[81] + F[58]*P[82]) + x82*F[58] + Q[13];
   Pnew[82] = x82;
   Pnew[83] = (x14*x14)*P[83] + Q[14];
   Pnew[84] = P[84] + Q[15];
   Pnew[85] = P[85] + Q[16];
   Pnew[86] = P[86] + Q[17];
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
   Pnew[12] = -x5*(P[1]*P[1]) + P[12];
   Pnew[13] = -x6*P[2] + P[13];
   Pnew[14] = -x6*P[3] + P[14];
   Pnew[15] = -x6*P[4] + P[15];
   Pnew[16] = -x6*P[5] + P[16];
   Pnew[17] = -x6*P[6] + P[17];
   Pnew[18] = -x6*P[7] + P[18];
   Pnew[19] = -x6*P[8] + P[19];
   Pnew[20] = -x6*P[9] + P[20];
   Pnew[21] = -x6*P[10] + P[21];
   Pnew[22] = -x6*P[11] + P[22];
   Pnew[23] = -x5*(P[2]*P[2]) + P[23];
   Pnew[24] = -x7*P[3] + P[24];
   Pnew[25] = -x7*P[4] + P[25];
   Pnew[26] = -x7*P[5] + P[26];
   Pnew[27] = -x7*P[6] + P[27];
   Pnew[28] = -x7*P[7] + P[28];
   Pnew[29] = -x7*P[8] + P[29];
   Pnew[30] = -x7*P[9] + P[30];
   Pnew[31] = -x7*P[10] + P[31];
   Pnew[32] = -x7*P[11] + P[32];
   Pnew[33] = -x5*(P[3]*P[3]) + P[33];
   Pnew[34] = -x8*P[4] + P[34];
   Pnew[35] = -x8*P[5] + P[35];
   Pnew[36] = -x8*P[6] + P[36];
   Pnew[37] = -x8*P[7] + P[37];
   Pnew[38] = -x8*P[8] + P[38];
   Pnew[39] = -x8*P[9] + P[39];
   Pnew[40] = -x8*P[10] + P[40];
   Pnew[41] = -x8*P[11] + P[41];
   Pnew[42] = -x5*(P[4]*P[4]) + P[42];
   Pnew[43] = -x9*P[5] + P[43];
   Pnew[44] = -x9*P[6] + P[44];
   Pnew[45] = -x9*P[7] + P[45];
   Pnew[46] = -x9*P[8] + P[46];
   Pnew[47] = -x9*P[9] + P[47];
   Pnew[48] = -x9*P[10] + P[48];
   Pnew[49] = -x9*P[11] + P[49];
   Pnew[50] = -x5*(P[5]*P[5]) + P[50];
   Pnew[51] = -x10*P[6] + P[51];
   Pnew[52] = -x10*P[7] + P[52];
   Pnew[53] = -x10*P[8] + P[53];
   Pnew[54] = -x10*P[9] + P[54];
   Pnew[55] = -x10*P[10] + P[55];
   Pnew[56] = -x10*P[11] + P[56];
   Pnew[57] = -x5*(P[6]*P[6]) + P[57];
   Pnew[58] = -x11*P[7] + P[58];
   Pnew[59] = -x11*P[8] + P[59];
   Pnew[60] = -x11*P[9] + P[60];
   Pnew[61] = -x11*P[10] + P[61];
   Pnew[62] = -x11*P[11] + P[62];
   Pnew[63] = -x5*(P[7]*P[7]) + P[63];
   Pnew[64] = -x12*P[8] + P[64];
   Pnew[65] = -x12*P[9] + P[65];
   Pnew[66] = -x12*P[10] + P[66];
   Pnew[67] = -x12*P[11] + P[67];
   Pnew[68] = -x5*(P[8]*P[8]) + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = -x5*(P[9]*P[9]) + P[72];
   Pnew[73] = P[73];
   Pnew[74] = -x5*(P[10]*P[10]) + P[74];
   Pnew[75] = P[75];
   Pnew[76] = P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = P[80];
   Pnew[81] = P[81];
   Pnew[82] = P[82];
   Pnew[83] = -x5*(P[11]*P[11]) + P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
}

void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, const float * restrict param,
        float *restrict xnew, float *restrict Pnew)
{
   const float x0 = 2*param[10];
   const float x1 = x[4]*x[7];
   const float x2 = x[5]*x[6];
   const float x3 = (x[4]*x[4]) - (x[7]*x[7]);
   const float x4 = (x[6]*x[6]);
   const float x5 = (x[5]*x[5]);
   const float x6 = x0*(x[4]*x[5] + x[6]*x[7]) - 2*(x1 - x2)*param[8] + (x3 + x4 - x5)*param[9] - mag[1];
   const float x7 = H[19]*P[4] + H[20]*P[5] + H[21]*P[6] + H[22]*P[7];
   const float x8 = H[19]*P[42] + H[20]*P[43] + H[21]*P[44] + H[22]*P[45];
   const float x9 = x8*H[23];
   const float x10 = H[19]*P[43] + H[20]*P[50] + H[21]*P[51] + H[22]*P[52];
   const float x11 = x10*H[24];
   const float x12 = H[19]*P[44] + H[20]*P[51] + H[21]*P[57] + H[22]*P[58];
   const float x13 = x12*H[25];
   const float x14 = H[19]*P[45] + H[20]*P[52] + H[21]*P[58] + H[22]*P[63];
   const float x15 = x14*H[26];
   const float x16 = H[23]*P[42] + H[24]*P[43] + H[25]*P[44] + H[26]*P[45];
   const float x17 = x16*H[19];
   const float x18 = H[23]*P[43] + H[24]*P[50] + H[25]*P[51] + H[26]*P[52];
   const float x19 = x18*H[20];
   const float x20 = H[23]*P[44] + H[24]*P[51] + H[25]*P[57] + H[26]*P[58];
   const float x21 = x20*H[21];
   const float x22 = H[23]*P[45] + H[24]*P[52] + H[25]*P[58] + H[26]*P[63];
   const float x23 = x22*H[22];
   const float x24 = x10*H[20] + x12*H[21] + x14*H[22] + x8*H[19] + R[7];
   const float x25 = x16*H[23] + x18*H[24] + x20*H[25] + x22*H[26] + R[8];
   const float x26 = 1.0F/(x24*x25 - (x11 + x13 + x15 + x9)*(x17 + x19 + x21 + x23));
   const float x27 = x26*(-x11 - x13 - x15 - x9);
   const float x28 = x27*x7;
   const float x29 = H[23]*P[4] + H[24]*P[5] + H[25]*P[6] + H[26]*P[7];
   const float x30 = x24*x26;
   const float x31 = x29*x30;
   const float x32 = -x0*(x[4]*x[6] - x[5]*x[7]) + 2*(x1 + x2)*param[9] + (x3 - x4 + x5)*param[8] - mag[0];
   const float x33 = x26*(-x17 - x19 - x21 - x23);
   const float x34 = x29*x33;
   const float x35 = x25*x26;
   const float x36 = x35*x7;
   const float x37 = H[19]*P[15] + H[20]*P[16] + H[21]*P[17] + H[22]*P[18];
   const float x38 = x27*x37;
   const float x39 = H[23]*P[15] + H[24]*P[16] + H[25]*P[17] + H[26]*P[18];
   const float x40 = x30*x39;
   const float x41 = x33*x39;
   const float x42 = x35*x37;
   const float x43 = H[19]*P[25] + H[20]*P[26] + H[21]*P[27] + H[22]*P[28];
   const float x44 = x27*x43;
   const float x45 = H[23]*P[25] + H[24]*P[26] + H[25]*P[27] + H[26]*P[28];
   const float x46 = x30*x45;
   const float x47 = x33*x45;
   const float x48 = x35*x43;
   const float x49 = H[19]*P[34] + H[20]*P[35] + H[21]*P[36] + H[22]*P[37];
   const float x50 = x27*x49;
   const float x51 = H[23]*P[34] + H[24]*P[35] + H[25]*P[36] + H[26]*P[37];
   const float x52 = x30*x51;
   const float x53 = x33*x51;
   const float x54 = x35*x49;
   const float x55 = x27*x8;
   const float x56 = x16*x30;
   const float x57 = x16*x33;
   const float x58 = x35*x8;
   const float x59 = x10*x27;
   const float x60 = x18*x30;
   const float x61 = x18*x33;
   const float x62 = x10*x35;
   const float x63 = x12*x27;
   const float x64 = x20*x30;
   const float x65 = x20*x33;
   const float x66 = x12*x35;
   const float x67 = x14*x27;
   const float x68 = x22*x30;
   const float x69 = x22*x33;
   const float x70 = x14*x35;
   const float x71 = H[19]*P[46] + H[20]*P[53] + H[21]*P[59] + H[22]*P[64];
   const float x72 = x27*x71;
   const float x73 = H[23]*P[46] + H[24]*P[53] + H[25]*P[59] + H[26]*P[64];
   const float x74 = x30*x73;
   const float x75 = x33*x73;
   const float x76 = x35*x71;
   const float x77 = H[19]*P[47] + H[20]*P[54] + H[21]*P[60] + H[22]*P[65];
   const float x78 = x27*x77;
   const float x79 = H[23]*P[47] + H[24]*P[54] + H[25]*P[60] + H[26]*P[65];
   const float x80 = x30*x79;
   const float x81 = x33*x79;
   const float x82 = x35*x77;
   const float x83 = H[19]*P[48] + H[20]*P[55] + H[21]*P[61] + H[22]*P[66];
   const float x84 = x27*x83;
   const float x85 = H[23]*P[48] + H[24]*P[55] + H[25]*P[61] + H[26]*P[66];
   const float x86 = x30*x85;
   const float x87 = x33*x85;
   const float x88 = x35*x83;
   const float x89 = H[19]*P[49] + H[20]*P[56] + H[21]*P[62] + H[22]*P[67];
   const float x90 = x27*x89;
   const float x91 = H[23]*P[49] + H[24]*P[56] + H[25]*P[62] + H[26]*P[67];
   const float x92 = x30*x91;
   const float x93 = x33*x91;
   const float x94 = x35*x89;
   const float x95 = x28 + x31;
   const float x96 = x34 + x36;
   const float x97 = -x95*H[23] - x96*H[19];
   const float x98 = -x95*H[24] - x96*H[20];
   const float x99 = -x95*H[25] - x96*H[21];
   const float x100 = -x95*H[26] - x96*H[22];
   const float x101 = x38 + x40;
   const float x102 = x41 + x42;
   const float x103 = -x101*H[23] - x102*H[19];
   const float x104 = -x101*H[24] - x102*H[20];
   const float x105 = -x101*H[25] - x102*H[21];
   const float x106 = -x101*H[26] - x102*H[22];
   const float x107 = x44 + x46;
   const float x108 = x47 + x48;
   const float x109 = -x107*H[23] - x108*H[19];
   const float x110 = -x107*H[24] - x108*H[20];
   const float x111 = -x107*H[25] - x108*H[21];
   const float x112 = -x107*H[26] - x108*H[22];
   const float x113 = x50 + x52;
   const float x114 = x53 + x54;
   const float x115 = -x113*H[23] - x114*H[19];
   const float x116 = -x113*H[24] - x114*H[20];
   const float x117 = -x113*H[25] - x114*H[21];
   const float x118 = -x113*H[26] - x114*H[22];
   const float x119 = x55 + x56;
   const float x120 = x57 + x58;
   const float x121 = -x119*H[24] - x120*H[20];
   const float x122 = -x119*H[25] - x120*H[21];
   const float x123 = -x119*H[26] - x120*H[22];
   const float x124 = -x119*H[23] - x120*H[19] + 1;
   const float x125 = x59 + x60;
   const float x126 = x61 + x62;
   const float x127 = -x125*H[23] - x126*H[19];
   const float x128 = -x125*H[25] - x126*H[21];
   const float x129 = -x125*H[26] - x126*H[22];
   const float x130 = -x125*H[24] - x126*H[20] + 1;
   const float x131 = x63 + x64;
   const float x132 = x65 + x66;
   const float x133 = -x131*H[23] - x132*H[19];
   const float x134 = -x131*H[24] - x132*H[20];
   const float x135 = -x131*H[26] - x132*H[22];
   const float x136 = -x131*H[25] - x132*H[21] + 1;
   const float x137 = x67 + x68;
   const float x138 = x69 + x70;
   const float x139 = -x137*H[23] - x138*H[19];
   const float x140 = -x137*H[24] - x138*H[20];
   const float x141 = -x137*H[25] - x138*H[21];
   const float x142 = -x137*H[26] - x138*H[22] + 1;
   const float x143 = x72 + x74;
   const float x144 = x75 + x76;
   const float x145 = x78 + x80;
   const float x146 = x81 + x82;
   const float x147 = x84 + x86;
   const float x148 = x87 + x88;
   const float x149 = x90 + x92;
   const float x150 = x93 + x94;
   xnew[0] = x32*(-x34 - x36) + x6*(-x28 - x31) + x[0];
   xnew[1] = x32*(-x41 - x42) + x6*(-x38 - x40) + x[1];
   xnew[2] = x32*(-x47 - x48) + x6*(-x44 - x46) + x[2];
   xnew[3] = x32*(-x53 - x54) + x6*(-x50 - x52) + x[3];
   xnew[4] = x32*(-x57 - x58) + x6*(-x55 - x56) + x[4];
   xnew[5] = x32*(-x61 - x62) + x6*(-x59 - x60) + x[5];
   xnew[6] = x32*(-x65 - x66) + x6*(-x63 - x64) + x[6];
   xnew[7] = x32*(-x69 - x70) + x6*(-x67 - x68) + x[7];
   xnew[8] = x32*(-x75 - x76) + x6*(-x72 - x74) + x[8];
   xnew[9] = x32*(-x81 - x82) + x6*(-x78 - x80) + x[9];
   xnew[10] = x32*(-x87 - x88) + x6*(-x84 - x86) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x32*(-x93 - x94) + x6*(-x90 - x92) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   Pnew[0] = x100*P[7] + x97*P[4] + x98*P[5] + x99*P[6] + P[0];
   Pnew[1] = x100*P[18] + x97*P[15] + x98*P[16] + x99*P[17] + P[1];
   Pnew[2] = x100*P[28] + x97*P[25] + x98*P[26] + x99*P[27] + P[2];
   Pnew[3] = x100*P[37] + x97*P[34] + x98*P[35] + x99*P[36] + P[3];
   Pnew[4] = x100*P[45] + x97*P[42] + x98*P[43] + x99*P[44] + P[4];
   Pnew[5] = x100*P[52] + x97*P[43] + x98*P[50] + x99*P[51] + P[5];
   Pnew[6] = x100*P[58] + x97*P[44] + x98*P[51] + x99*P[57] + P[6];
   Pnew[7] = x100*P[63] + x97*P[45] + x98*P[52] + x99*P[58] + P[7];
   Pnew[8] = x100*P[64] + x97*P[46] + x98*P[53] + x99*P[59] + P[8];
   Pnew[9] = x100*P[65] + x97*P[47] + x98*P[54] + x99*P[60] + P[9];
   Pnew[10] = x100*P[66] + x97*P[48] + x98*P[55] + x99*P[61] + P[10];
   Pnew[11] = x100*P[67] + x97*P[49] + x98*P[56] + x99*P[62] + P[11];
   Pnew[12] = x103*P[15] + x104*P[16] + x105*P[17] + x106*P[18] + P[12];
   Pnew[13] = x103*P[25] + x104*P[26] + x105*P[27] + x106*P[28] + P[13];
   Pnew[14] = x103*P[34] + x104*P[35] + x105*P[36] + x106*P[37] + P[14];
   Pnew[15] = x103*P[42] + x104*P[43] + x105*P[44] + x106*P[45] + P[15];
   Pnew[16] = x103*P[43] + x104*P[50] + x105*P[51] + x106*P[52] + P[16];
   Pnew[17] = x103*P[44] + x104*P[51] + x105*P[57] + x106*P[58] + P[17];
   Pnew[18] = x103*P[45] + x104*P[52] + x105*P[58] + x106*P[63] + P[18];
   Pnew[19] = x103*P[46] + x104*P[53] + x105*P[59] + x106*P[64] + P[19];
   Pnew[20] = x103*P[47] + x104*P[54] + x105*P[60] + x106*P[65] + P[20];
   Pnew[21] = x103*P[48] + x104*P[55] + x105*P[61] + x106*P[66] + P[21];
   Pnew[22] = x103*P[49] + x104*P[56] + x105*P[62] + x106*P[67] + P[22];
   Pnew[23] = x109*P[25] + x110*P[26] + x111*P[27] + x112*P[28] + P[23];
   Pnew[24] = x109*P[34] + x110*P[35] + x111*P[36] + x112*P[37] + P[24];
   Pnew[25] = x109*P[42] + x110*P[43] + x111*P[44] + x112*P[45] + P[25];
   Pnew[26] = x109*P[43] + x110*P[50] + x111*P[51] + x112*P[52] + P[26];
   Pnew[27] = x109*P[44] + x110*P[51] + x111*P[57] + x112*P[58] + P[27];
   Pnew[28] = x109*P[45] + x110*P[52] + x111*P[58] + x112*P[63] + P[28];
   Pnew[29] = x109*P[46] + x110*P[53] + x111*P[59] + x112*P[64] + P[29];
   Pnew[30] = x109*P[47] + x110*P[54] + x111*P[60] + x112*P[65] + P[30];
   Pnew[31] = x109*P[48] + x110*P[55] + x111*P[61] + x112*P[66] + P[31];
   Pnew[32] = x109*P[49] + x110*P[56] + x111*P[62] + x112*P[67] + P[32];
   Pnew[33] = x115*P[34] + x116*P[35] + x117*P[36] + x118*P[37] + P[33];
   Pnew[34] = x115*P[42] + x116*P[43] + x117*P[44] + x118*P[45] + P[34];
   Pnew[35] = x115*P[43] + x116*P[50] + x117*P[51] + x118*P[52] + P[35];
   Pnew[36] = x115*P[44] + x116*P[51] + x117*P[57] + x118*P[58] + P[36];
   Pnew[37] = x115*P[45] + x116*P[52] + x117*P[58] + x118*P[63] + P[37];
   Pnew[38] = x115*P[46] + x116*P[53] + x117*P[59] + x118*P[64] + P[38];
   Pnew[39] = x115*P[47] + x116*P[54] + x117*P[60] + x118*P[65] + P[39];
   Pnew[40] = x115*P[48] + x116*P[55] + x117*P[61] + x118*P[66] + P[40];
   Pnew[41] = x115*P[49] + x116*P[56] + x117*P[62] + x118*P[67] + P[41];
   Pnew[42] = x121*P[43] + x122*P[44] + x123*P[45] + x124*P[42];
   Pnew[43] = x121*P[50] + x122*P[51] + x123*P[52] + x124*P[43];
   Pnew[44] = x121*P[51] + x122*P[57] + x123*P[58] + x124*P[44];
   Pnew[45] = x121*P[52] + x122*P[58] + x123*P[63] + x124*P[45];
   Pnew[46] = x121*P[53] + x122*P[59] + x123*P[64] + x124*P[46];
   Pnew[47] = x121*P[54] + x122*P[60] + x123*P[65] + x124*P[47];
   Pnew[48] = x121*P[55] + x122*P[61] + x123*P[66] + x124*P[48];
   Pnew[49] = x121*P[56] + x122*P[62] + x123*P[67] + x124*P[49];
   Pnew[50] = x127*P[43] + x128*P[51] + x129*P[52] + x130*P[50];
   Pnew[51] = x127*P[44] + x128*P[57] + x129*P[58] + x130*P[51];
   Pnew[52] = x127*P[45] + x128*P[58] + x129*P[63] + x130*P[52];
   Pnew[53] = x127*P[46] + x128*P[59] + x129*P[64] + x130*P[53];
   Pnew[54] = x127*P[47] + x128*P[60] + x129*P[65] + x130*P[54];
   Pnew[55] = x127*P[48] + x128*P[61] + x129*P[66] + x130*P[55];
   Pnew[56] = x127*P[49] + x128*P[62] + x129*P[67] + x130*P[56];
   Pnew[57] = x133*P[44] + x134*P[51] + x135*P[58] + x136*P[57];
   Pnew[58] = x133*P[45] + x134*P[52] + x135*P[63] + x136*P[58];
   Pnew[59] = x133*P[46] + x134*P[53] + x135*P[64] + x136*P[59];
   Pnew[60] = x133*P[47] + x134*P[54] + x135*P[65] + x136*P[60];
   Pnew[61] = x133*P[48] + x134*P[55] + x135*P[66] + x136*P[61];
   Pnew[62] = x133*P[49] + x134*P[56] + x135*P[67] + x136*P[62];
   Pnew[63] = x139*P[45] + x140*P[52] + x141*P[58] + x142*P[63];
   Pnew[64] = x139*P[46] + x140*P[53] + x141*P[59] + x142*P[64];
   Pnew[65] = x139*P[47] + x140*P[54] + x141*P[60] + x142*P[65];
   Pnew[66] = x139*P[48] + x140*P[55] + x141*P[61] + x142*P[66];
   Pnew[67] = x139*P[49] + x140*P[56] + x141*P[62] + x142*P[67];
   Pnew[68] = (-x143*H[23] - x144*H[19])*P[46] + (-x143*H[24] - x144*H[20])*P[53] + (-x143*H[25] - x144*H[21])*P[59] + (-x143*H[26] - x144*H[22])*P[64] + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = (-x145*H[23] - x146*H[19])*P[47] + (-x145*H[24] - x146*H[20])*P[54] + (-x145*H[25] - x146*H[21])*P[60] + (-x145*H[26] - x146*H[22])*P[65] + P[72];
   Pnew[73] = P[73];
   Pnew[74] = (-x147*H[23] - x148*H[19])*P[48] + (-x147*H[24] - x148*H[20])*P[55] + (-x147*H[25] - x148*H[21])*P[61] + (-x147*H[26] - x148*H[22])*P[66] + P[74];
   Pnew[75] = P[75];
   Pnew[76] = P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = P[80];
   Pnew[81] = P[81];
   Pnew[82] = P[82];
   Pnew[83] = (-x149*H[23] - x150*H[19])*P[49] + (-x149*H[24] - x150*H[20])*P[56] + (-x149*H[25] - x150*H[21])*P[62] + (-x149*H[26] - x150*H[22])*P[67] + P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
}

void gyro_correction(const float * restrict x, const float * restrict P, const float * restrict gyro,
	const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = (H[16]*H[16]);
   const float x1 = x0*P[68];
   const float x2 = 1.0F/(x1 + R[4]);
   const float x3 = x2*(-gyro[0] + x[8])*H[16];
   const float x4 = (H[17]*H[17]);
   const float x5 = x4*P[72];
   const float x6 = 1.0F/(x5 + R[5]);
   const float x7 = x6*(-gyro[1] + x[9])*H[17];
   const float x8 = (H[18]*H[18]);
   const float x9 = x8*P[74];
   const float x10 = 1.0F/(x9 + R[6]);
   const float x11 = x10*(-gyro[2] + x[10])*H[18];
   const float x12 = x0*x2;
   const float x13 = x4*x6;
   const float x14 = x10*x8;
   const float x15 = x0*x2*P[8];
   const float x16 = x4*x6*P[9];
   const float x17 = x10*x8*P[10];
   const float x18 = x1*x2;
   const float x19 = x5*x6;
   const float x20 = x10*x9;
   const float x21 = x0*x2*P[19];
   const float x22 = x4*x6*P[20];
   const float x23 = x10*x8*P[21];
   const float x24 = x0*x2*P[29];
   const float x25 = x4*x6*P[30];
   const float x26 = x10*x8*P[31];
   const float x27 = x0*x2*P[38];
   const float x28 = x4*x6*P[39];
   const float x29 = x10*x8*P[40];
   const float x30 = x0*x2*P[46];
   const float x31 = x4*x6*P[47];
   const float x32 = x10*x8*P[48];
   const float x33 = x0*x2*P[53];
   const float x34 = x4*x6*P[54];
   const float x35 = x10*x8*P[55];
   const float x36 = -x18 + 1;
   const float x37 = -x19 + 1;
   const float x38 = -x20 + 1;
   xnew[0] = -x11*P[10] - x3*P[8] - x7*P[9] + x[0];
   xnew[1] = -x11*P[21] - x3*P[19] - x7*P[20] + x[1];
   xnew[2] = -x11*P[31] - x3*P[29] - x7*P[30] + x[2];
   xnew[3] = -x11*P[40] - x3*P[38] - x7*P[39] + x[3];
   xnew[4] = -x11*P[48] - x3*P[46] - x7*P[47] + x[4];
   xnew[5] = -x11*P[55] - x3*P[53] - x7*P[54] + x[5];
   xnew[6] = -x11*P[61] - x3*P[59] - x7*P[60] + x[6];
   xnew[7] = -x11*P[66] - x3*P[64] - x7*P[65] + x[7];
   xnew[8] = -x3*P[68] + x[8];
   xnew[9] = -x7*P[72] + x[9];
   xnew[10] = -x11*P[74] + x[10];
   xnew[11] = -x3*P[69] + x[11];
   xnew[12] = -x3*P[70] + x[12];
   xnew[13] = -x11*P[75] + x[13];
   xnew[14] = x[14];
   xnew[15] = -x3*P[71] + x[15];
   xnew[16] = -x7*P[73] + x[16];
   xnew[17] = -x11*P[76] + x[17];
   Pnew[0] = -x12*(P[8]*P[8]) - x13*(P[9]*P[9]) - x14*(P[10]*P[10]) + P[0];
   Pnew[1] = -x15*P[19] - x16*P[20] - x17*P[21] + P[1];
   Pnew[2] = -x15*P[29] - x16*P[30] - x17*P[31] + P[2];
   Pnew[3] = -x15*P[38] - x16*P[39] - x17*P[40] + P[3];
   Pnew[4] = -x15*P[46] - x16*P[47] - x17*P[48] + P[4];
   Pnew[5] = -x15*P[53] - x16*P[54] - x17*P[55] + P[5];
   Pnew[6] = -x15*P[59] - x16*P[60] - x17*P[61] + P[6];
   Pnew[7] = -x15*P[64] - x16*P[65] - x17*P[66] + P[7];
   Pnew[8] = -x18*P[8] + P[8];
   Pnew[9] = -x19*P[9] + P[9];
   Pnew[10] = -x20*P[10] + P[10];
   Pnew[11] = P[11];
   Pnew[12] = -x12*(P[19]*P[19]) - x13*(P[20]*P[20]) - x14*(P[21]*P[21]) + P[12];
   Pnew[13] = -x21*P[29] - x22*P[30] - x23*P[31] + P[13];
   Pnew[14] = -x21*P[38] - x22*P[39] - x23*P[40] + P[14];
   Pnew[15] = -x21*P[46] - x22*P[47] - x23*P[48] + P[15];
   Pnew[16] = -x21*P[53] - x22*P[54] - x23*P[55] + P[16];
   Pnew[17] = -x21*P[59] - x22*P[60] - x23*P[61] + P[17];
   Pnew[18] = -x21*P[64] - x22*P[65] - x23*P[66] + P[18];
   Pnew[19] = -x18*P[19] + P[19];
   Pnew[20] = -x19*P[20] + P[20];
   Pnew[21] = -x20*P[21] + P[21];
   Pnew[22] = P[22];
   Pnew[23] = -x12*(P[29]*P[29]) - x13*(P[30]*P[30]) - x14*(P[31]*P[31]) + P[23];
   Pnew[24] = -x24*P[38] - x25*P[39] - x26*P[40] + P[24];
   Pnew[25] = -x24*P[46] - x25*P[47] - x26*P[48] + P[25];
   Pnew[26] = -x24*P[53] - x25*P[54] - x26*P[55] + P[26];
   Pnew[27] = -x24*P[59] - x25*P[60] - x26*P[61] + P[27];
   Pnew[28] = -x24*P[64] - x25*P[65] - x26*P[66] + P[28];
   Pnew[29] = -x18*P[29] + P[29];
   Pnew[30] = -x19*P[30] + P[30];
   Pnew[31] = -x20*P[31] + P[31];
   Pnew[32] = P[32];
   Pnew[33] = -x12*(P[38]*P[38]) - x13*(P[39]*P[39]) - x14*(P[40]*P[40]) + P[33];
   Pnew[34] = -x27*P[46] - x28*P[47] - x29*P[48] + P[34];
   Pnew[35] = -x27*P[53] - x28*P[54] - x29*P[55] + P[35];
   Pnew[36] = -x27*P[59] - x28*P[60] - x29*P[61] + P[36];
   Pnew[37] = -x27*P[64] - x28*P[65] - x29*P[66] + P[37];
   Pnew[38] = -x18*P[38] + P[38];
   Pnew[39] = -x19*P[39] + P[39];
   Pnew[40] = -x20*P[40] + P[40];
   Pnew[41] = P[41];
   Pnew[42] = -x12*(P[46]*P[46]) - x13*(P[47]*P[47]) - x14*(P[48]*P[48]) + P[42];
   Pnew[43] = -x30*P[53] - x31*P[54] - x32*P[55] + P[43];
   Pnew[44] = -x30*P[59] - x31*P[60] - x32*P[61] + P[44];
   Pnew[45] = -x30*P[64] - x31*P[65] - x32*P[66] + P[45];
   Pnew[46] = -x18*P[46] + P[46];
   Pnew[47] = -x19*P[47] + P[47];
   Pnew[48] = -x20*P[48] + P[48];
   Pnew[49] = P[49];
   Pnew[50] = -x12*(P[53]*P[53]) - x13*(P[54]*P[54]) - x14*(P[55]*P[55]) + P[50];
   Pnew[51] = -x33*P[59] - x34*P[60] - x35*P[61] + P[51];
   Pnew[52] = -x33*P[64] - x34*P[65] - x35*P[66] + P[52];
   Pnew[53] = -x18*P[53] + P[53];
   Pnew[54] = -x19*P[54] + P[54];
   Pnew[55] = -x20*P[55] + P[55];
   Pnew[56] = P[56];
   Pnew[57] = -x12*(P[59]*P[59]) - x13*(P[60]*P[60]) - x14*(P[61]*P[61]) + P[57];
   Pnew[58] = -x12*P[59]*P[64] - x13*P[60]*P[65] - x14*P[61]*P[66] + P[58];
   Pnew[59] = -x18*P[59] + P[59];
   Pnew[60] = -x19*P[60] + P[60];
   Pnew[61] = -x20*P[61] + P[61];
   Pnew[62] = P[62];
   Pnew[63] = -x12*(P[64]*P[64]) - x13*(P[65]*P[65]) - x14*(P[66]*P[66]) + P[63];
   Pnew[64] = -x18*P[64] + P[64];
   Pnew[65] = -x19*P[65] + P[65];
   Pnew[66] = -x20*P[66] + P[66];
   Pnew[67] = P[67];
   Pnew[68] = x36*P[68];
   Pnew[69] = x36*P[69];
   Pnew[70] = x36*P[70];
   Pnew[71] = x36*P[71];
   Pnew[72] = x37*P[72];
   Pnew[73] = x37*P[73];
   Pnew[74] = x38*P[74];
   Pnew[75] = x38*P[75];
   Pnew[76] = x38*P[76];
   Pnew[77] = -x12*(P[69]*P[69]) + P[77];
   Pnew[78] = -x12*P[69]*P[71] + P[78];
   Pnew[79] = -x12*(P[70]*P[70]) + P[79];
   Pnew[80] = P[80];
   Pnew[81] = -x14*(P[75]*P[75]) + P[81];
   Pnew[82] = -x14*P[75]*P[76] + P[82];
   Pnew[83] = P[83];
   Pnew[84] = -x12*(P[71]*P[71]) + P[84];
   Pnew[85] = -x13*(P[73]*P[73]) + P[85];
   Pnew[86] = -x14*(P[76]*P[76]) + P[86];
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
   const float x6 = -(x0*(x[4]*x[5] + x[6]*x[7]) - 2*(x1 - x2)*x[1] + (x3 + x4 - x5)*x[2])*param[7] - accel[1];
   const float x7 = H[1]*P[1] + H[2]*P[2] + H[3]*P[3] + H[4]*P[4] + H[5]*P[5] + H[6]*P[6] + H[7]*P[7];
   const float x8 = H[1]*P[22];
   const float x9 = H[2]*P[32];
   const float x10 = H[3]*P[41];
   const float x11 = H[4]*P[49];
   const float x12 = H[5]*P[56];
   const float x13 = H[6]*P[62];
   const float x14 = H[7]*P[67];
   const float x15 = x10 + x11 + x12 + x13 + x14 + x8 + x9;
   const float x16 = x15*H[15];
   const float x17 = H[8]*P[22];
   const float x18 = H[9]*P[32];
   const float x19 = H[10]*P[41];
   const float x20 = H[11]*P[49];
   const float x21 = H[12]*P[56];
   const float x22 = H[13]*P[62];
   const float x23 = H[14]*P[67];
   const float x24 = x17*H[15] + x18*H[15] + x19*H[15] + x20*H[15] + x21*H[15] + x22*H[15] + x23*H[15];
   const float x25 = x16*x24;
   const float x26 = (H[15]*H[15])*P[83] + R[3];
   const float x27 = H[1]*P[12] + H[2]*P[13] + H[3]*P[14] + H[4]*P[15] + H[5]*P[16] + H[6]*P[17] + H[7]*P[18];
   const float x28 = H[1]*P[13] + H[2]*P[23] + H[3]*P[24] + H[4]*P[25] + H[5]*P[26] + H[6]*P[27] + H[7]*P[28];
   const float x29 = H[1]*P[14] + H[2]*P[24] + H[3]*P[33] + H[4]*P[34] + H[5]*P[35] + H[6]*P[36] + H[7]*P[37];
   const float x30 = H[1]*P[15] + H[2]*P[25] + H[3]*P[34] + H[4]*P[42] + H[5]*P[43] + H[6]*P[44] + H[7]*P[45];
   const float x31 = H[1]*P[16] + H[2]*P[26] + H[3]*P[35] + H[4]*P[43] + H[5]*P[50] + H[6]*P[51] + H[7]*P[52];
   const float x32 = H[1]*P[17] + H[2]*P[27] + H[3]*P[36] + H[4]*P[44] + H[5]*P[51] + H[6]*P[57] + H[7]*P[58];
   const float x33 = H[1]*P[18] + H[2]*P[28] + H[3]*P[37] + H[4]*P[45] + H[5]*P[52] + H[6]*P[58] + H[7]*P[63];
   const float x34 = x27*H[8] + x28*H[9] + x29*H[10] + x30*H[11] + x31*H[12] + x32*H[13] + x33*H[14];
   const float x35 = x26*x34;
   const float x36 = x17 + x18 + x19 + x20 + x21 + x22 + x23;
   const float x37 = x36*H[15];
   const float x38 = x10*H[15] + x11*H[15] + x12*H[15] + x13*H[15] + x14*H[15] + x8*H[15] + x9*H[15];
   const float x39 = x37*x38;
   const float x40 = H[8]*P[12] + H[9]*P[13] + H[10]*P[14] + H[11]*P[15] + H[12]*P[16] + H[13]*P[17] + H[14]*P[18];
   const float x41 = H[8]*P[13] + H[9]*P[23] + H[10]*P[24] + H[11]*P[25] + H[12]*P[26] + H[13]*P[27] + H[14]*P[28];
   const float x42 = H[8]*P[14] + H[9]*P[24] + H[10]*P[33] + H[11]*P[34] + H[12]*P[35] + H[13]*P[36] + H[14]*P[37];
   const float x43 = H[8]*P[15] + H[9]*P[25] + H[10]*P[34] + H[11]*P[42] + H[12]*P[43] + H[13]*P[44] + H[14]*P[45];
   const float x44 = H[8]*P[16] + H[9]*P[26] + H[10]*P[35] + H[11]*P[43] + H[12]*P[50] + H[13]*P[51] + H[14]*P[52];
   const float x45 = H[8]*P[17] + H[9]*P[27] + H[10]*P[36] + H[11]*P[44] + H[12]*P[51] + H[13]*P[57] + H[14]*P[58];
   const float x46 = H[8]*P[18] + H[9]*P[28] + H[10]*P[37] + H[11]*P[45] + H[12]*P[52] + H[13]*P[58] + H[14]*P[63];
   const float x47 = x40*H[1] + x41*H[2] + x42*H[3] + x43*H[4] + x44*H[5] + x45*H[6] + x46*H[7];
   const float x48 = x16*x38;
   const float x49 = x40*H[8] + x41*H[9] + x42*H[10] + x43*H[11] + x44*H[12] + x45*H[13] + x46*H[14] + R[2];
   const float x50 = x24*x37;
   const float x51 = x27*H[1] + x28*H[2] + x29*H[3] + x30*H[4] + x31*H[5] + x32*H[6] + x33*H[7] + R[1];
   const float x52 = x26*x51;
   const float x53 = 1.0F/(x25*x47 + x34*x39 - x35*x47 - x48*x49 + x49*x52 - x50*x51);
   const float x54 = x53*(x25 - x35);
   const float x55 = x54*x7;
   const float x56 = H[8]*P[1] + H[9]*P[2] + H[10]*P[3] + H[11]*P[4] + H[12]*P[5] + H[13]*P[6] + H[14]*P[7];
   const float x57 = x53*(-x48 + x52);
   const float x58 = x56*x57;
   const float x59 = x53*(-x24*x51 + x34*x38)*H[15];
   const float x60 = x59*P[11];
   const float x61 = -(-x0*(x[4]*x[6] - x[5]*x[7]) + 2*(x1 + x2)*x[2] + (x3 - x4 + x5)*x[1])*param[7] - accel[0];
   const float x62 = x53*(-x26*x47 + x39);
   const float x63 = x56*x62;
   const float x64 = x53*(x26*x49 - x50);
   const float x65 = x64*x7;
   const float x66 = x53*(x24*x47 - x38*x49)*H[15];
   const float x67 = x66*P[11];
   const float x68 = -accel[2] - param[0]*param[5]*x[14];
   const float x69 = x53*(-x16*x49 + x34*x37);
   const float x70 = x69*x7;
   const float x71 = x53*(x16*x47 - x37*x51);
   const float x72 = x56*x71;
   const float x73 = x53*(-x34*x47 + x49*x51)*H[15];
   const float x74 = x73*P[11];
   const float x75 = x27*x54;
   const float x76 = x40*x57;
   const float x77 = x59*P[22];
   const float x78 = x40*x62;
   const float x79 = x27*x64;
   const float x80 = x66*P[22];
   const float x81 = x27*x69;
   const float x82 = x40*x71;
   const float x83 = x73*P[22];
   const float x84 = x28*x54;
   const float x85 = x41*x57;
   const float x86 = x59*P[32];
   const float x87 = x41*x62;
   const float x88 = x28*x64;
   const float x89 = x66*P[32];
   const float x90 = x28*x69;
   const float x91 = x41*x71;
   const float x92 = x73*P[32];
   const float x93 = x29*x54;
   const float x94 = x42*x57;
   const float x95 = x59*P[41];
   const float x96 = x42*x62;
   const float x97 = x29*x64;
   const float x98 = x66*P[41];
   const float x99 = x29*x69;
   const float x100 = x42*x71;
   const float x101 = x73*P[41];
   const float x102 = x30*x54;
   const float x103 = x43*x57;
   const float x104 = x59*P[49];
   const float x105 = x43*x62;
   const float x106 = x30*x64;
   const float x107 = x66*P[49];
   const float x108 = x30*x69;
   const float x109 = x43*x71;
   const float x110 = x73*P[49];
   const float x111 = x31*x54;
   const float x112 = x44*x57;
   const float x113 = x59*P[56];
   const float x114 = x44*x62;
   const float x115 = x31*x64;
   const float x116 = x66*P[56];
   const float x117 = x31*x69;
   const float x118 = x44*x71;
   const float x119 = x73*P[56];
   const float x120 = x32*x54;
   const float x121 = x45*x57;
   const float x122 = x59*P[62];
   const float x123 = x45*x62;
   const float x124 = x32*x64;
   const float x125 = x66*P[62];
   const float x126 = x32*x69;
   const float x127 = x45*x71;
   const float x128 = x73*P[62];
   const float x129 = x33*x54;
   const float x130 = x46*x57;
   const float x131 = x59*P[67];
   const float x132 = x46*x62;
   const float x133 = x33*x64;
   const float x134 = x66*P[67];
   const float x135 = x33*x69;
   const float x136 = x46*x71;
   const float x137 = x73*P[67];
   const float x138 = H[1]*P[19] + H[2]*P[29] + H[3]*P[38] + H[4]*P[46] + H[5]*P[53] + H[6]*P[59] + H[7]*P[64];
   const float x139 = x138*x54;
   const float x140 = H[8]*P[19] + H[9]*P[29] + H[10]*P[38] + H[11]*P[46] + H[12]*P[53] + H[13]*P[59] + H[14]*P[64];
   const float x141 = x140*x57;
   const float x142 = x140*x62;
   const float x143 = x138*x64;
   const float x144 = H[1]*P[20] + H[2]*P[30] + H[3]*P[39] + H[4]*P[47] + H[5]*P[54] + H[6]*P[60] + H[7]*P[65];
   const float x145 = x144*x54;
   const float x146 = H[8]*P[20] + H[9]*P[30] + H[10]*P[39] + H[11]*P[47] + H[12]*P[54] + H[13]*P[60] + H[14]*P[65];
   const float x147 = x146*x57;
   const float x148 = x146*x62;
   const float x149 = x144*x64;
   const float x150 = H[1]*P[21] + H[2]*P[31] + H[3]*P[40] + H[4]*P[48] + H[5]*P[55] + H[6]*P[61] + H[7]*P[66];
   const float x151 = x150*x54;
   const float x152 = H[8]*P[21] + H[9]*P[31] + H[10]*P[40] + H[11]*P[48] + H[12]*P[55] + H[13]*P[61] + H[14]*P[66];
   const float x153 = x152*x57;
   const float x154 = x152*x62;
   const float x155 = x150*x64;
   const float x156 = x15*x54;
   const float x157 = x36*x57;
   const float x158 = x59*P[83];
   const float x159 = x36*x62;
   const float x160 = x15*x64;
   const float x161 = x66*P[83];
   const float x162 = x15*x69;
   const float x163 = x36*x71;
   const float x164 = x73*P[83];
   const float x165 = (x70 + x72 + x74)*H[15];
   const float x166 = x55 + x58 + x60;
   const float x167 = x63 + x65 + x67;
   const float x168 = -x166*H[8] - x167*H[1];
   const float x169 = -x166*H[9] - x167*H[2];
   const float x170 = -x166*H[10] - x167*H[3];
   const float x171 = -x166*H[11] - x167*H[4];
   const float x172 = -x166*H[12] - x167*H[5];
   const float x173 = -x166*H[13] - x167*H[6];
   const float x174 = -x166*H[14] - x167*H[7];
   const float x175 = (x81 + x82 + x83)*H[15];
   const float x176 = x75 + x76 + x77;
   const float x177 = x78 + x79 + x80;
   const float x178 = -x176*H[9] - x177*H[2];
   const float x179 = -x176*H[10] - x177*H[3];
   const float x180 = -x176*H[11] - x177*H[4];
   const float x181 = -x176*H[12] - x177*H[5];
   const float x182 = -x176*H[13] - x177*H[6];
   const float x183 = -x176*H[14] - x177*H[7];
   const float x184 = -x176*H[8] - x177*H[1] + 1;
   const float x185 = (x90 + x91 + x92)*H[15];
   const float x186 = x84 + x85 + x86;
   const float x187 = x87 + x88 + x89;
   const float x188 = -x186*H[8] - x187*H[1];
   const float x189 = -x186*H[10] - x187*H[3];
   const float x190 = -x186*H[11] - x187*H[4];
   const float x191 = -x186*H[12] - x187*H[5];
   const float x192 = -x186*H[13] - x187*H[6];
   const float x193 = -x186*H[14] - x187*H[7];
   const float x194 = -x186*H[9] - x187*H[2] + 1;
   const float x195 = (x100 + x101 + x99)*H[15];
   const float x196 = x93 + x94 + x95;
   const float x197 = x96 + x97 + x98;
   const float x198 = -x196*H[8] - x197*H[1];
   const float x199 = -x196*H[9] - x197*H[2];
   const float x200 = -x196*H[11] - x197*H[4];
   const float x201 = -x196*H[12] - x197*H[5];
   const float x202 = -x196*H[13] - x197*H[6];
   const float x203 = -x196*H[14] - x197*H[7];
   const float x204 = -x196*H[10] - x197*H[3] + 1;
   const float x205 = (x108 + x109 + x110)*H[15];
   const float x206 = x102 + x103 + x104;
   const float x207 = x105 + x106 + x107;
   const float x208 = -x206*H[8] - x207*H[1];
   const float x209 = -x206*H[9] - x207*H[2];
   const float x210 = -x206*H[10] - x207*H[3];
   const float x211 = -x206*H[12] - x207*H[5];
   const float x212 = -x206*H[13] - x207*H[6];
   const float x213 = -x206*H[14] - x207*H[7];
   const float x214 = -x206*H[11] - x207*H[4] + 1;
   const float x215 = (x117 + x118 + x119)*H[15];
   const float x216 = x111 + x112 + x113;
   const float x217 = x114 + x115 + x116;
   const float x218 = -x216*H[8] - x217*H[1];
   const float x219 = -x216*H[9] - x217*H[2];
   const float x220 = -x216*H[10] - x217*H[3];
   const float x221 = -x216*H[11] - x217*H[4];
   const float x222 = -x216*H[13] - x217*H[6];
   const float x223 = -x216*H[14] - x217*H[7];
   const float x224 = -x216*H[12] - x217*H[5] + 1;
   const float x225 = (x126 + x127 + x128)*H[15];
   const float x226 = x120 + x121 + x122;
   const float x227 = x123 + x124 + x125;
   const float x228 = -x226*H[8] - x227*H[1];
   const float x229 = -x226*H[9] - x227*H[2];
   const float x230 = -x226*H[10] - x227*H[3];
   const float x231 = -x226*H[11] - x227*H[4];
   const float x232 = -x226*H[12] - x227*H[5];
   const float x233 = -x226*H[14] - x227*H[7];
   const float x234 = -x226*H[13] - x227*H[6] + 1;
   const float x235 = (x135 + x136 + x137)*H[15];
   const float x236 = x129 + x130 + x131;
   const float x237 = x132 + x133 + x134;
   const float x238 = -x236*H[8] - x237*H[1];
   const float x239 = -x236*H[9] - x237*H[2];
   const float x240 = -x236*H[10] - x237*H[3];
   const float x241 = -x236*H[11] - x237*H[4];
   const float x242 = -x236*H[12] - x237*H[5];
   const float x243 = -x236*H[13] - x237*H[6];
   const float x244 = -x236*H[14] - x237*H[7] + 1;
   const float x245 = x139 + x141;
   const float x246 = x142 + x143;
   const float x247 = x145 + x147;
   const float x248 = x148 + x149;
   const float x249 = x151 + x153;
   const float x250 = x154 + x155;
   const float x251 = x156 + x157 + x158;
   const float x252 = x159 + x160 + x161;
   xnew[0] = x6*(-x55 - x58 - x60) + x61*(-x63 - x65 - x67) + x68*(-x70 - x72 - x74) + x[0];
   xnew[1] = x6*(-x75 - x76 - x77) + x61*(-x78 - x79 - x80) + x68*(-x81 - x82 - x83) + x[1];
   xnew[2] = x6*(-x84 - x85 - x86) + x61*(-x87 - x88 - x89) + x68*(-x90 - x91 - x92) + x[2];
   xnew[3] = x6*(-x93 - x94 - x95) + x61*(-x96 - x97 - x98) + x68*(-x100 - x101 - x99) + x[3];
   xnew[4] = x6*(-x102 - x103 - x104) + x61*(-x105 - x106 - x107) + x68*(-x108 - x109 - x110) + x[4];
   xnew[5] = x6*(-x111 - x112 - x113) + x61*(-x114 - x115 - x116) + x68*(-x117 - x118 - x119) + x[5];
   xnew[6] = x6*(-x120 - x121 - x122) + x61*(-x123 - x124 - x125) + x68*(-x126 - x127 - x128) + x[6];
   xnew[7] = x6*(-x129 - x130 - x131) + x61*(-x132 - x133 - x134) + x68*(-x135 - x136 - x137) + x[7];
   xnew[8] = x6*(-x139 - x141) + x61*(-x142 - x143) + x68*(-x138*x69 - x140*x71) + x[8];
   xnew[9] = x6*(-x145 - x147) + x61*(-x148 - x149) + x68*(-x144*x69 - x146*x71) + x[9];
   xnew[10] = x6*(-x151 - x153) + x61*(-x154 - x155) + x68*(-x150*x69 - x152*x71) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x6*(-x156 - x157 - x158) + x61*(-x159 - x160 - x161) + x68*(-x162 - x163 - x164) + x[14];
   xnew[15] = x[15];
   xnew[16] = x[16];
   xnew[17] = x[17];
   Pnew[0] = -x165*P[11] + x168*P[1] + x169*P[2] + x170*P[3] + x171*P[4] + x172*P[5] + x173*P[6] + x174*P[7] + P[0];
   Pnew[1] = -x165*P[22] + x168*P[12] + x169*P[13] + x170*P[14] + x171*P[15] + x172*P[16] + x173*P[17] + x174*P[18] + P[1];
   Pnew[2] = -x165*P[32] + x168*P[13] + x169*P[23] + x170*P[24] + x171*P[25] + x172*P[26] + x173*P[27] + x174*P[28] + P[2];
   Pnew[3] = -x165*P[41] + x168*P[14] + x169*P[24] + x170*P[33] + x171*P[34] + x172*P[35] + x173*P[36] + x174*P[37] + P[3];
   Pnew[4] = -x165*P[49] + x168*P[15] + x169*P[25] + x170*P[34] + x171*P[42] + x172*P[43] + x173*P[44] + x174*P[45] + P[4];
   Pnew[5] = -x165*P[56] + x168*P[16] + x169*P[26] + x170*P[35] + x171*P[43] + x172*P[50] + x173*P[51] + x174*P[52] + P[5];
   Pnew[6] = -x165*P[62] + x168*P[17] + x169*P[27] + x170*P[36] + x171*P[44] + x172*P[51] + x173*P[57] + x174*P[58] + P[6];
   Pnew[7] = -x165*P[67] + x168*P[18] + x169*P[28] + x170*P[37] + x171*P[45] + x172*P[52] + x173*P[58] + x174*P[63] + P[7];
   Pnew[8] = x168*P[19] + x169*P[29] + x170*P[38] + x171*P[46] + x172*P[53] + x173*P[59] + x174*P[64] + P[8];
   Pnew[9] = x168*P[20] + x169*P[30] + x170*P[39] + x171*P[47] + x172*P[54] + x173*P[60] + x174*P[65] + P[9];
   Pnew[10] = x168*P[21] + x169*P[31] + x170*P[40] + x171*P[48] + x172*P[55] + x173*P[61] + x174*P[66] + P[10];
   Pnew[11] = -x165*P[83] + x168*P[22] + x169*P[32] + x170*P[41] + x171*P[49] + x172*P[56] + x173*P[62] + x174*P[67] + P[11];
   Pnew[12] = -x175*P[22] + x178*P[13] + x179*P[14] + x180*P[15] + x181*P[16] + x182*P[17] + x183*P[18] + x184*P[12];
   Pnew[13] = -x175*P[32] + x178*P[23] + x179*P[24] + x180*P[25] + x181*P[26] + x182*P[27] + x183*P[28] + x184*P[13];
   Pnew[14] = -x175*P[41] + x178*P[24] + x179*P[33] + x180*P[34] + x181*P[35] + x182*P[36] + x183*P[37] + x184*P[14];
   Pnew[15] = -x175*P[49] + x178*P[25] + x179*P[34] + x180*P[42] + x181*P[43] + x182*P[44] + x183*P[45] + x184*P[15];
   Pnew[16] = -x175*P[56] + x178*P[26] + x179*P[35] + x180*P[43] + x181*P[50] + x182*P[51] + x183*P[52] + x184*P[16];
   Pnew[17] = -x175*P[62] + x178*P[27] + x179*P[36] + x180*P[44] + x181*P[51] + x182*P[57] + x183*P[58] + x184*P[17];
   Pnew[18] = -x175*P[67] + x178*P[28] + x179*P[37] + x180*P[45] + x181*P[52] + x182*P[58] + x183*P[63] + x184*P[18];
   Pnew[19] = x178*P[29] + x179*P[38] + x180*P[46] + x181*P[53] + x182*P[59] + x183*P[64] + x184*P[19];
   Pnew[20] = x178*P[30] + x179*P[39] + x180*P[47] + x181*P[54] + x182*P[60] + x183*P[65] + x184*P[20];
   Pnew[21] = x178*P[31] + x179*P[40] + x180*P[48] + x181*P[55] + x182*P[61] + x183*P[66] + x184*P[21];
   Pnew[22] = -x175*P[83] + x178*P[32] + x179*P[41] + x180*P[49] + x181*P[56] + x182*P[62] + x183*P[67] + x184*P[22];
   Pnew[23] = -x185*P[32] + x188*P[13] + x189*P[24] + x190*P[25] + x191*P[26] + x192*P[27] + x193*P[28] + x194*P[23];
   Pnew[24] = -x185*P[41] + x188*P[14] + x189*P[33] + x190*P[34] + x191*P[35] + x192*P[36] + x193*P[37] + x194*P[24];
   Pnew[25] = -x185*P[49] + x188*P[15] + x189*P[34] + x190*P[42] + x191*P[43] + x192*P[44] + x193*P[45] + x194*P[25];
   Pnew[26] = -x185*P[56] + x188*P[16] + x189*P[35] + x190*P[43] + x191*P[50] + x192*P[51] + x193*P[52] + x194*P[26];
   Pnew[27] = -x185*P[62] + x188*P[17] + x189*P[36] + x190*P[44] + x191*P[51] + x192*P[57] + x193*P[58] + x194*P[27];
   Pnew[28] = -x185*P[67] + x188*P[18] + x189*P[37] + x190*P[45] + x191*P[52] + x192*P[58] + x193*P[63] + x194*P[28];
   Pnew[29] = x188*P[19] + x189*P[38] + x190*P[46] + x191*P[53] + x192*P[59] + x193*P[64] + x194*P[29];
   Pnew[30] = x188*P[20] + x189*P[39] + x190*P[47] + x191*P[54] + x192*P[60] + x193*P[65] + x194*P[30];
   Pnew[31] = x188*P[21] + x189*P[40] + x190*P[48] + x191*P[55] + x192*P[61] + x193*P[66] + x194*P[31];
   Pnew[32] = -x185*P[83] + x188*P[22] + x189*P[41] + x190*P[49] + x191*P[56] + x192*P[62] + x193*P[67] + x194*P[32];
   Pnew[33] = -x195*P[41] + x198*P[14] + x199*P[24] + x200*P[34] + x201*P[35] + x202*P[36] + x203*P[37] + x204*P[33];
   Pnew[34] = -x195*P[49] + x198*P[15] + x199*P[25] + x200*P[42] + x201*P[43] + x202*P[44] + x203*P[45] + x204*P[34];
   Pnew[35] = -x195*P[56] + x198*P[16] + x199*P[26] + x200*P[43] + x201*P[50] + x202*P[51] + x203*P[52] + x204*P[35];
   Pnew[36] = -x195*P[62] + x198*P[17] + x199*P[27] + x200*P[44] + x201*P[51] + x202*P[57] + x203*P[58] + x204*P[36];
   Pnew[37] = -x195*P[67] + x198*P[18] + x199*P[28] + x200*P[45] + x201*P[52] + x202*P[58] + x203*P[63] + x204*P[37];
   Pnew[38] = x198*P[19] + x199*P[29] + x200*P[46] + x201*P[53] + x202*P[59] + x203*P[64] + x204*P[38];
   Pnew[39] = x198*P[20] + x199*P[30] + x200*P[47] + x201*P[54] + x202*P[60] + x203*P[65] + x204*P[39];
   Pnew[40] = x198*P[21] + x199*P[31] + x200*P[48] + x201*P[55] + x202*P[61] + x203*P[66] + x204*P[40];
   Pnew[41] = -x195*P[83] + x198*P[22] + x199*P[32] + x200*P[49] + x201*P[56] + x202*P[62] + x203*P[67] + x204*P[41];
   Pnew[42] = -x205*P[49] + x208*P[15] + x209*P[25] + x210*P[34] + x211*P[43] + x212*P[44] + x213*P[45] + x214*P[42];
   Pnew[43] = -x205*P[56] + x208*P[16] + x209*P[26] + x210*P[35] + x211*P[50] + x212*P[51] + x213*P[52] + x214*P[43];
   Pnew[44] = -x205*P[62] + x208*P[17] + x209*P[27] + x210*P[36] + x211*P[51] + x212*P[57] + x213*P[58] + x214*P[44];
   Pnew[45] = -x205*P[67] + x208*P[18] + x209*P[28] + x210*P[37] + x211*P[52] + x212*P[58] + x213*P[63] + x214*P[45];
   Pnew[46] = x208*P[19] + x209*P[29] + x210*P[38] + x211*P[53] + x212*P[59] + x213*P[64] + x214*P[46];
   Pnew[47] = x208*P[20] + x209*P[30] + x210*P[39] + x211*P[54] + x212*P[60] + x213*P[65] + x214*P[47];
   Pnew[48] = x208*P[21] + x209*P[31] + x210*P[40] + x211*P[55] + x212*P[61] + x213*P[66] + x214*P[48];
   Pnew[49] = -x205*P[83] + x208*P[22] + x209*P[32] + x210*P[41] + x211*P[56] + x212*P[62] + x213*P[67] + x214*P[49];
   Pnew[50] = -x215*P[56] + x218*P[16] + x219*P[26] + x220*P[35] + x221*P[43] + x222*P[51] + x223*P[52] + x224*P[50];
   Pnew[51] = -x215*P[62] + x218*P[17] + x219*P[27] + x220*P[36] + x221*P[44] + x222*P[57] + x223*P[58] + x224*P[51];
   Pnew[52] = -x215*P[67] + x218*P[18] + x219*P[28] + x220*P[37] + x221*P[45] + x222*P[58] + x223*P[63] + x224*P[52];
   Pnew[53] = x218*P[19] + x219*P[29] + x220*P[38] + x221*P[46] + x222*P[59] + x223*P[64] + x224*P[53];
   Pnew[54] = x218*P[20] + x219*P[30] + x220*P[39] + x221*P[47] + x222*P[60] + x223*P[65] + x224*P[54];
   Pnew[55] = x218*P[21] + x219*P[31] + x220*P[40] + x221*P[48] + x222*P[61] + x223*P[66] + x224*P[55];
   Pnew[56] = -x215*P[83] + x218*P[22] + x219*P[32] + x220*P[41] + x221*P[49] + x222*P[62] + x223*P[67] + x224*P[56];
   Pnew[57] = -x225*P[62] + x228*P[17] + x229*P[27] + x230*P[36] + x231*P[44] + x232*P[51] + x233*P[58] + x234*P[57];
   Pnew[58] = -x225*P[67] + x228*P[18] + x229*P[28] + x230*P[37] + x231*P[45] + x232*P[52] + x233*P[63] + x234*P[58];
   Pnew[59] = x228*P[19] + x229*P[29] + x230*P[38] + x231*P[46] + x232*P[53] + x233*P[64] + x234*P[59];
   Pnew[60] = x228*P[20] + x229*P[30] + x230*P[39] + x231*P[47] + x232*P[54] + x233*P[65] + x234*P[60];
   Pnew[61] = x228*P[21] + x229*P[31] + x230*P[40] + x231*P[48] + x232*P[55] + x233*P[66] + x234*P[61];
   Pnew[62] = -x225*P[83] + x228*P[22] + x229*P[32] + x230*P[41] + x231*P[49] + x232*P[56] + x233*P[67] + x234*P[62];
   Pnew[63] = -x235*P[67] + x238*P[18] + x239*P[28] + x240*P[37] + x241*P[45] + x242*P[52] + x243*P[58] + x244*P[63];
   Pnew[64] = x238*P[19] + x239*P[29] + x240*P[38] + x241*P[46] + x242*P[53] + x243*P[59] + x244*P[64];
   Pnew[65] = x238*P[20] + x239*P[30] + x240*P[39] + x241*P[47] + x242*P[54] + x243*P[60] + x244*P[65];
   Pnew[66] = x238*P[21] + x239*P[31] + x240*P[40] + x241*P[48] + x242*P[55] + x243*P[61] + x244*P[66];
   Pnew[67] = -x235*P[83] + x238*P[22] + x239*P[32] + x240*P[41] + x241*P[49] + x242*P[56] + x243*P[62] + x244*P[67];
   Pnew[68] = (-x245*H[8] - x246*H[1])*P[19] + (-x245*H[9] - x246*H[2])*P[29] + (-x245*H[10] - x246*H[3])*P[38] + (-x245*H[11] - x246*H[4])*P[46] + (-x245*H[12] - x246*H[5])*P[53] + (-x245*H[13] - x246*H[6])*P[59] + (-x245*H[14] - x246*H[7])*P[64] + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = (-x247*H[8] - x248*H[1])*P[20] + (-x247*H[9] - x248*H[2])*P[30] + (-x247*H[10] - x248*H[3])*P[39] + (-x247*H[11] - x248*H[4])*P[47] + (-x247*H[12] - x248*H[5])*P[54] + (-x247*H[13] - x248*H[6])*P[60] + (-x247*H[14] - x248*H[7])*P[65] + P[72];
   Pnew[73] = P[73];
   Pnew[74] = (-x249*H[8] - x250*H[1])*P[21] + (-x249*H[9] - x250*H[2])*P[31] + (-x249*H[10] - x250*H[3])*P[40] + (-x249*H[11] - x250*H[4])*P[48] + (-x249*H[12] - x250*H[5])*P[55] + (-x249*H[13] - x250*H[6])*P[61] + (-x249*H[14] - x250*H[7])*P[66] + P[74];
   Pnew[75] = P[75];
   Pnew[76] = P[76];
   Pnew[77] = P[77];
   Pnew[78] = P[78];
   Pnew[79] = P[79];
   Pnew[80] = P[80];
   Pnew[81] = P[81];
   Pnew[82] = P[82];
   Pnew[83] = (-x251*H[8] - x252*H[1])*P[22] + (-x251*H[9] - x252*H[2])*P[32] + (-x251*H[10] - x252*H[3])*P[41] + (-x251*H[11] - x252*H[4])*P[49] + (-x251*H[12] - x252*H[5])*P[56] + (-x251*H[13] - x252*H[6])*P[62] + (-x251*H[14] - x252*H[7])*P[67] + (-(x162 + x163 + x164)*H[15] + 1)*P[83];
   Pnew[84] = P[84];
   Pnew[85] = P[85];
   Pnew[86] = P[86];
}

void covariance_init(float *Pnew) {
	float P0[NUMX] = {1.0e7f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		              1.0e-9f, 1.0e-9f, 1.0e-9f};
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
   Pnew[12] = P0[1];
   Pnew[13] = 0;
   Pnew[14] = 0;
   Pnew[15] = 0;
   Pnew[16] = 0;
   Pnew[17] = 0;
   Pnew[18] = 0;
   Pnew[19] = 0;
   Pnew[20] = 0;
   Pnew[21] = 0;
   Pnew[22] = 0;
   Pnew[23] = P0[2];
   Pnew[24] = 0;
   Pnew[25] = 0;
   Pnew[26] = 0;
   Pnew[27] = 0;
   Pnew[28] = 0;
   Pnew[29] = 0;
   Pnew[30] = 0;
   Pnew[31] = 0;
   Pnew[32] = 0;
   Pnew[33] = P0[3];
   Pnew[34] = 0;
   Pnew[35] = 0;
   Pnew[36] = 0;
   Pnew[37] = 0;
   Pnew[38] = 0;
   Pnew[39] = 0;
   Pnew[40] = 0;
   Pnew[41] = 0;
   Pnew[42] = P0[4];
   Pnew[43] = 0;
   Pnew[44] = 0;
   Pnew[45] = 0;
   Pnew[46] = 0;
   Pnew[47] = 0;
   Pnew[48] = 0;
   Pnew[49] = 0;
   Pnew[50] = P0[5];
   Pnew[51] = 0;
   Pnew[52] = 0;
   Pnew[53] = 0;
   Pnew[54] = 0;
   Pnew[55] = 0;
   Pnew[56] = 0;
   Pnew[57] = P0[6];
   Pnew[58] = 0;
   Pnew[59] = 0;
   Pnew[60] = 0;
   Pnew[61] = 0;
   Pnew[62] = 0;
   Pnew[63] = P0[7];
   Pnew[64] = 0;
   Pnew[65] = 0;
   Pnew[66] = 0;
   Pnew[67] = 0;
   Pnew[68] = P0[8];
   Pnew[69] = 0;
   Pnew[70] = 0;
   Pnew[71] = 0;
   Pnew[72] = P0[9];
   Pnew[73] = 0;
   Pnew[74] = P0[10];
   Pnew[75] = 0;
   Pnew[76] = 0;
   Pnew[77] = P0[11];
   Pnew[78] = 0;
   Pnew[79] = P0[12];
   Pnew[80] = 0;
   Pnew[81] = P0[13];
   Pnew[82] = 0;
   Pnew[83] = P0[14];
   Pnew[84] = P0[15];
   Pnew[85] = P0[16];
   Pnew[86] = P0[17];
}

/**
 * @}
 * @}
 */
