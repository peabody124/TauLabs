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
#define NUMX 15  // number of state variables
#define NUMF 56  // number of non-zero variables stored in the F matrix
#define NUMH 27  // number of non-zero variables stored in the H matrix
#define NUMP 74  // number of non-zero variables in the upper triangular of covariance
#define NUMO 9   // number of total outputs

// Private methods
void covariance_init(float *Pnew);
void state_prediction(const float * restrict x, const float *u, float Ts, const float *restrict param, float * restrict xnew);
void covariance_prediction(const float *restrict P, float Ts, const float * restrict F, const float * restrict Q, float * restrict Pnew);
void linearize_FH(const float * restrict x, const float * restrict param, float * restrict F, float * restrict H);
void baro_correction(const float *restrict x, const float *restrict P, float baro,
        const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew);
void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew);
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
		float beta1;
		float tau;
		float mu;
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
	for (uint32_t i = 0; i < NUMX; i++)
		qcins_state->Q[i] = 1e-4;
	qcins_state->Q[NUMX-1] = 1;
	qcins_state->Q[NUMX-2] = 1;
	qcins_state->Q[NUMX-3] = 1;
	qcins_state->Q[NUMX-4] = 1;

	// Store the observation noises
	qcins_state->R[0] = 100;
	qcins_state->R[1] = qcins_state->R[2] = qcins_state->R[3] = 1e6;
	qcins_state->R[4] = qcins_state->R[5] = qcins_state->R[6] = 1000;
	qcins_state->R[7] = qcins_state->R[8] = 200;
	
	return true;
}

bool qcins_alloc(uintptr_t *qcins_handle)
{
	struct qcins_state *qcins_state = (struct qcins_state *) PIOS_malloc(sizeof(struct qcins_state));
	if (qcins_state == NULL)
		return false;

	qcins_state->magic = QCINS_STATE_MAGIC;

	qcins_state->params.g = 9.81f;
	qcins_state->params.beta1 = 10000.0f / 180.0f * 3.1415f;
	qcins_state->params.tau = 0.050f;
	qcins_state->params.mu = 1;

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


bool qcins_set_gains(uintptr_t qcins_handle, const float gains_new[4])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	qcins_state->params.beta1 = expf(gains_new[0]) * DEG2RAD;
	// TODO: expand to have beta for each axis
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
	linearize_FH(qcins_state->x, &qcins_state->params.g, qcins_state->F, qcins_state->H);
	state_prediction(qcins_state->x, u, Ts, &qcins_state->params.g, qcins_state->xnew);
	covariance_prediction(qcins_state->P, Ts, qcins_state->F, qcins_state->Q, qcins_state->Pnew);
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

	mag_correction(qcins_state->x, qcins_state->P, mag_norm, qcins_state->H, qcins_state->R, qcins_state->xnew, qcins_state->Pnew);
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
   const float x0 = 2*param[0]*x[14];
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
   const float x17 = x16*param[3];
   const float x18 = -x12;
   const float x19 = x10 + x14 + x18 + x9;
   const float x20 = x3 + 2*x[5]*x[6];
   const float x21 = x1 - x2;
   const float x22 = x19*x[1] + x20*x[2] - x21*x5;
   const float x23 = x22*param[3];
   const float x24 = 2*param[3];
   const float x25 = (1.0F/2.0F)*x[5];
   const float x26 = (1.0F/2.0F)*x[6];
   const float x27 = (1.0F/2.0F)*x[7];
   const float x28 = (1.0F/2.0F)*x[4];
   const float x29 = Ts*param[1];
   const float x30 = Ts/param[2];
   xnew[0] = Ts*x[3] + x[0];
   xnew[1] = Ts*(-x0*(x1 + x2) + x17*x4 - x19*x23) + x[1];
   xnew[2] = Ts*(x0*(x6 - x7) - x15*x17 - x20*x23) + x[2];
   xnew[3] = Ts*(-x16*x24*x8 + x21*x22*x24 - (x11 + x13 + x18)*param[0]*x[14] + param[0]) + x[3];
   xnew[4] = Ts*(-x25*x[8] - x26*x[9] - x27*x[10]) + x[4];
   xnew[5] = Ts*(x26*x[10] - x27*x[9] + x28*x[8]) + x[5];
   xnew[6] = Ts*(-x25*x[10] + x27*x[8] + x28*x[9]) + x[6];
   xnew[7] = Ts*(x25*x[9] - x26*x[8] + x28*x[10]) + x[7];
   xnew[8] = x29*x[11] + x[8];
   xnew[9] = x29*x[12] + x[9];
   xnew[10] = x29*x[13] + x[10];
   xnew[11] = x30*(u[0] - x[11]) + x[11];
   xnew[12] = x30*(u[1] - x[12]) + x[12];
   xnew[13] = x30*(u[2] - x[13]) + x[13];
   xnew[14] = x30*(u[3] - x[14]) + x[14];
}

void linearize_FH(const float * restrict x, const float * restrict param, float * restrict F, float * restrict H)
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
   const float x12 = x4*param[3];
   const float x13 = x5 - x6;
   const float x14 = x10 + x13 + x7;
   const float x15 = x1 + x3;
   const float x16 = x15*param[3];
   const float x17 = -x11*x16 - x12*x14;
   const float x18 = x0*x[5];
   const float x19 = 2*x[6];
   const float x20 = x19*x[7];
   const float x21 = x18 + x20;
   const float x22 = x21*param[3];
   const float x23 = x0*x[6];
   const float x24 = x2*x[7];
   const float x25 = -x23 + x24;
   const float x26 = x25*param[3];
   const float x27 = -x11*x26 - x22*x4;
   const float x28 = 2*x[7];
   const float x29 = x0*x[2] + x2*x[3] - x28*x[1];
   const float x30 = x0*x[1];
   const float x31 = x28*x[2];
   const float x32 = x19*x[3];
   const float x33 = x30 + x31 - x32;
   const float x34 = x11*param[3];
   const float x35 = param[0]*x[14];
   const float x36 = x19*x35;
   const float x37 = (x14*x[2] + x21*x[3] + x4*x[1])*param[3];
   const float x38 = x28*x37;
   const float x39 = (x11*x[1] + x15*x[2] + x25*x[3])*param[3];
   const float x40 = x0*x39;
   const float x41 = -x36 + x38 - x40;
   const float x42 = x19*x[1];
   const float x43 = x2*x[2];
   const float x44 = x0*x[3];
   const float x45 = x42 - x43 + x44;
   const float x46 = x19*x[2] + x2*x[1] + x28*x[3];
   const float x47 = -x19*x37 - x2*x39 - x28*x35;
   const float x48 = -x42 + x43 - x44;
   const float x49 = x0*x35;
   const float x50 = x2*x37;
   const float x51 = x19*x39;
   const float x52 = -x49 - x50 + x51;
   const float x53 = x2*x35;
   const float x54 = -x30 - x31 + x32;
   const float x55 = x0*x37;
   const float x56 = x28*x39;
   const float x57 = -x14*x22 - x15*x26;
   const float x58 = x14*param[3];
   const float x59 = x53 - x55 - x56;
   const float x60 = (1.0F/2.0F)*x[8];
   const float x61 = -x60;
   const float x62 = (1.0F/2.0F)*x[9];
   const float x63 = -x62;
   const float x64 = (1.0F/2.0F)*x[10];
   const float x65 = -x64;
   const float x66 = (1.0F/2.0F)*x[5];
   const float x67 = -x66;
   const float x68 = (1.0F/2.0F)*x[6];
   const float x69 = -x68;
   const float x70 = (1.0F/2.0F)*x[7];
   const float x71 = -x70;
   const float x72 = (1.0F/2.0F)*x[4];
   const float x73 = -1/param[2];
   const float x74 = -x46*param[3];
   const float x75 = -x29*param[3];
   const float x76 = -x28;
   F[0] = 1;
   F[1] = -(x11*x11)*param[3] - (x4*x4)*param[3];
   F[2] = x17;
   F[3] = x27;
   F[4] = -x12*x29 - x33*x34 + x41;
   F[5] = -x12*x45 - x34*x46 + x47;
   F[6] = -x12*x46 - x34*x48 + x52;
   F[7] = -x12*x54 - x29*x34 - x53 + x55 + x56;
   F[8] = -(x23 + x24)*param[0];
   F[9] = x17;
   F[10] = -(x14*x14)*param[3] - (x15*x15)*param[3];
   F[11] = x57;
   F[12] = -x16*x33 - x29*x58 + x59;
   F[13] = -x16*x46 - x45*x58 + x49 + x50 - x51;
   F[14] = -x16*x48 - x46*x58 + x47;
   F[15] = -x16*x29 + x41 - x54*x58;
   F[16] = -(-x18 + x20)*param[0];
   F[17] = x27;
   F[18] = x57;
   F[19] = -(x21*x21)*param[3] - (x25*x25)*param[3];
   F[20] = -x22*x29 - x26*x33 + x52;
   F[21] = -x22*x45 - x26*x46 + x59;
   F[22] = -x22*x46 - x26*x48 + x36 - x38 + x40;
   F[23] = -x22*x54 - x26*x29 + x47;
   F[24] = -(x13 + x8 + x9)*param[0];
   F[25] = x61;
   F[26] = x63;
   F[27] = x65;
   F[28] = x67;
   F[29] = x69;
   F[30] = x71;
   F[31] = x60;
   F[32] = x64;
   F[33] = x63;
   F[34] = x72;
   F[35] = x71;
   F[36] = x68;
   F[37] = x62;
   F[38] = x65;
   F[39] = x60;
   F[40] = x70;
   F[41] = x72;
   F[42] = x67;
   F[43] = x64;
   F[44] = x62;
   F[45] = x61;
   F[46] = x69;
   F[47] = x66;
   F[48] = x72;
   F[49] = param[1];
   F[50] = param[1];
   F[51] = param[1];
   F[52] = x73;
   F[53] = x73;
   F[54] = x73;
   F[55] = x73;
   H[0] = 1;
   H[1] = -x34;
   H[2] = -x16;
   H[3] = -x26;
   H[4] = -x33*param[3];
   H[5] = x74;
   H[6] = -x48*param[3];
   H[7] = x75;
   H[8] = -x12;
   H[9] = -x58;
   H[10] = -x22;
   H[11] = x75;
   H[12] = -x45*param[3];
   H[13] = x74;
   H[14] = -x54*param[3];
   H[15] = -param[0];
   H[16] = 1;
   H[17] = 1;
   H[18] = 1;
   H[19] = x0;
   H[20] = x2;
   H[21] = -x19;
   H[22] = x76;
   H[23] = x76;
   H[24] = x19;
   H[25] = x2;
   H[26] = -x0;
}

// TODO: rework F definition so that Ts is already multiplied into all the terms
// to reduce the memory requirement (currently duplicates F)
void covariance_prediction(const float *restrict P, float Ts, const float * restrict F, const float * restrict Q, float * restrict Pnew)
{
   const float x0 = Ts*F[0];
   const float x1 = x0*P[33] + P[3];
   const float x2 = Ts*F[2];
   const float x3 = x0*P[24] + P[2];
   const float x4 = Ts*F[3];
   const float x5 = Ts*F[4];
   const float x6 = x0*P[34] + P[4];
   const float x7 = Ts*F[5];
   const float x8 = x0*P[35] + P[5];
   const float x9 = Ts*F[6];
   const float x10 = x0*P[36] + P[6];
   const float x11 = Ts*F[7];
   const float x12 = x0*P[37] + P[7];
   const float x13 = Ts*F[8];
   const float x14 = x0*P[41] + P[11];
   const float x15 = Ts*F[1] + 1;
   const float x16 = x0*P[14] + P[1];
   const float x17 = Ts*F[9];
   const float x18 = Ts*F[11];
   const float x19 = Ts*F[12];
   const float x20 = Ts*F[13];
   const float x21 = Ts*F[14];
   const float x22 = Ts*F[15];
   const float x23 = Ts*F[16];
   const float x24 = Ts*F[10] + 1;
   const float x25 = Ts*F[17];
   const float x26 = Ts*F[18];
   const float x27 = Ts*F[20];
   const float x28 = Ts*F[21];
   const float x29 = Ts*F[22];
   const float x30 = Ts*F[23];
   const float x31 = Ts*F[24];
   const float x32 = Ts*F[19] + 1;
   const float x33 = Ts*F[25];
   const float x34 = Ts*F[26];
   const float x35 = Ts*F[27];
   const float x36 = Ts*F[28];
   const float x37 = x0*P[38] + P[8];
   const float x38 = Ts*F[29];
   const float x39 = x0*P[39] + P[9];
   const float x40 = Ts*F[30];
   const float x41 = x0*P[40] + P[10];
   const float x42 = Ts*F[31];
   const float x43 = Ts*F[32];
   const float x44 = Ts*F[33];
   const float x45 = Ts*F[34];
   const float x46 = Ts*F[35];
   const float x47 = Ts*F[36];
   const float x48 = Ts*F[37];
   const float x49 = Ts*F[38];
   const float x50 = Ts*F[39];
   const float x51 = Ts*F[40];
   const float x52 = Ts*F[41];
   const float x53 = Ts*F[42];
   const float x54 = Ts*F[43];
   const float x55 = Ts*F[44];
   const float x56 = Ts*F[45];
   const float x57 = Ts*F[46];
   const float x58 = Ts*F[47];
   const float x59 = Ts*F[48];
   const float x60 = Ts*F[55] + 1;
   const float x61 = x13*P[73] + x15*P[22] + x2*P[32] + x4*P[41];
   const float x62 = x11*P[45] + x15*P[15] + x2*P[25] + x4*P[34] + x5*P[42] + x7*P[43] + x9*P[44];
   const float x63 = x11*P[51] + x15*P[16] + x2*P[26] + x4*P[35] + x5*P[43] + x7*P[49] + x9*P[50];
   const float x64 = x11*P[56] + x15*P[17] + x2*P[27] + x4*P[36] + x5*P[44] + x7*P[50] + x9*P[55];
   const float x65 = x11*P[60] + x15*P[18] + x2*P[28] + x4*P[37] + x5*P[45] + x7*P[51] + x9*P[56];
   const float x66 = x11*P[28] + x13*P[32] + x15*P[13] + x2*P[23] + x4*P[24] + x5*P[25] + x7*P[26] + x9*P[27];
   const float x67 = x11*P[37] + x13*P[41] + x15*P[14] + x2*P[24] + x4*P[33] + x5*P[34] + x7*P[35] + x9*P[36];
   const float x68 = x11*P[18] + x13*P[22] + x15*P[12] + x2*P[13] + x4*P[14] + x5*P[15] + x7*P[16] + x9*P[17];
   const float x69 = x11*P[61] + x15*P[19] + x2*P[29] + x4*P[38] + x5*P[46] + x7*P[52] + x9*P[57];
   const float x70 = x11*P[62] + x15*P[20] + x2*P[30] + x4*P[39] + x5*P[47] + x7*P[53] + x9*P[58];
   const float x71 = x11*P[63] + x15*P[21] + x2*P[31] + x4*P[40] + x5*P[48] + x7*P[54] + x9*P[59];
   const float x72 = x17*P[22] + x18*P[41] + x23*P[73] + x24*P[32];
   const float x73 = x17*P[15] + x18*P[34] + x19*P[42] + x20*P[43] + x21*P[44] + x22*P[45] + x24*P[25];
   const float x74 = x17*P[16] + x18*P[35] + x19*P[43] + x20*P[49] + x21*P[50] + x22*P[51] + x24*P[26];
   const float x75 = x17*P[17] + x18*P[36] + x19*P[44] + x20*P[50] + x21*P[55] + x22*P[56] + x24*P[27];
   const float x76 = x17*P[18] + x18*P[37] + x19*P[45] + x20*P[51] + x21*P[56] + x22*P[60] + x24*P[28];
   const float x77 = x17*P[12] + x18*P[14] + x19*P[15] + x20*P[16] + x21*P[17] + x22*P[18] + x23*P[22] + x24*P[13];
   const float x78 = x17*P[14] + x18*P[33] + x19*P[34] + x20*P[35] + x21*P[36] + x22*P[37] + x23*P[41] + x24*P[24];
   const float x79 = x17*P[13] + x18*P[24] + x19*P[25] + x20*P[26] + x21*P[27] + x22*P[28] + x23*P[32] + x24*P[23];
   const float x80 = x17*P[19] + x18*P[38] + x19*P[46] + x20*P[52] + x21*P[57] + x22*P[61] + x24*P[29];
   const float x81 = x17*P[20] + x18*P[39] + x19*P[47] + x20*P[53] + x21*P[58] + x22*P[62] + x24*P[30];
   const float x82 = x17*P[21] + x18*P[40] + x19*P[48] + x20*P[54] + x21*P[59] + x22*P[63] + x24*P[31];
   const float x83 = x25*P[22] + x26*P[32] + x31*P[73] + x32*P[41];
   const float x84 = x25*P[15] + x26*P[25] + x27*P[42] + x28*P[43] + x29*P[44] + x30*P[45] + x32*P[34];
   const float x85 = x25*P[16] + x26*P[26] + x27*P[43] + x28*P[49] + x29*P[50] + x30*P[51] + x32*P[35];
   const float x86 = x25*P[17] + x26*P[27] + x27*P[44] + x28*P[50] + x29*P[55] + x30*P[56] + x32*P[36];
   const float x87 = x25*P[18] + x26*P[28] + x27*P[45] + x28*P[51] + x29*P[56] + x30*P[60] + x32*P[37];
   const float x88 = x25*P[19] + x26*P[29] + x27*P[46] + x28*P[52] + x29*P[57] + x30*P[61] + x32*P[38];
   const float x89 = x25*P[20] + x26*P[30] + x27*P[47] + x28*P[53] + x29*P[58] + x30*P[62] + x32*P[39];
   const float x90 = x25*P[21] + x26*P[31] + x27*P[48] + x28*P[54] + x29*P[59] + x30*P[63] + x32*P[40];
   const float x91 = x33*P[43] + x34*P[44] + x35*P[45] + x36*P[46] + x38*P[47] + x40*P[48] + P[42];
   const float x92 = x33*P[52] + x34*P[57] + x35*P[61] + x36*P[64] + P[46];
   const float x93 = x33*P[53] + x34*P[58] + x35*P[62] + x38*P[66] + P[47];
   const float x94 = x33*P[54] + x34*P[59] + x35*P[63] + x40*P[68] + P[48];
   const float x95 = x33*P[49] + x34*P[50] + x35*P[51] + x36*P[52] + x38*P[53] + x40*P[54] + P[43];
   const float x96 = x33*P[50] + x34*P[55] + x35*P[56] + x36*P[57] + x38*P[58] + x40*P[59] + P[44];
   const float x97 = x33*P[51] + x34*P[56] + x35*P[60] + x36*P[61] + x38*P[62] + x40*P[63] + P[45];
   const float x98 = (Ts*Ts);
   const float x99 = x98*F[49]*P[65];
   const float x100 = x98*F[50]*P[67];
   const float x101 = x98*F[51]*P[69];
   const float x102 = x42*P[43] + x43*P[50] + x44*P[51] + x45*P[52] + x46*P[53] + x47*P[54] + P[49];
   const float x103 = x42*P[46] + x43*P[57] + x44*P[61] + x45*P[64] + P[52];
   const float x104 = x42*P[47] + x43*P[58] + x44*P[62] + x46*P[66] + P[53];
   const float x105 = x42*P[48] + x43*P[59] + x44*P[63] + x47*P[68] + P[54];
   const float x106 = x42*P[42] + x43*P[44] + x44*P[45] + x45*P[46] + x46*P[47] + x47*P[48] + P[43];
   const float x107 = x42*P[44] + x43*P[55] + x44*P[56] + x45*P[57] + x46*P[58] + x47*P[59] + P[50];
   const float x108 = x42*P[45] + x43*P[56] + x44*P[60] + x45*P[61] + x46*P[62] + x47*P[63] + P[51];
   const float x109 = x48*P[44] + x49*P[50] + x50*P[56] + x51*P[57] + x52*P[58] + x53*P[59] + P[55];
   const float x110 = x48*P[46] + x49*P[52] + x50*P[61] + x51*P[64] + P[57];
   const float x111 = x48*P[47] + x49*P[53] + x50*P[62] + x52*P[66] + P[58];
   const float x112 = x48*P[48] + x49*P[54] + x50*P[63] + x53*P[68] + P[59];
   const float x113 = x48*P[42] + x49*P[43] + x50*P[45] + x51*P[46] + x52*P[47] + x53*P[48] + P[44];
   const float x114 = x48*P[43] + x49*P[49] + x50*P[51] + x51*P[52] + x52*P[53] + x53*P[54] + P[50];
   const float x115 = x48*P[45] + x49*P[51] + x50*P[60] + x51*P[61] + x52*P[62] + x53*P[63] + P[56];
   const float x116 = x54*P[46] + x55*P[52] + x56*P[57] + x57*P[64] + P[61];
   const float x117 = x54*P[47] + x55*P[53] + x56*P[58] + x58*P[66] + P[62];
   const float x118 = x54*P[48] + x55*P[54] + x56*P[59] + x59*P[68] + P[63];
   const float x119 = Ts*F[49];
   const float x120 = x119*P[70] + P[65];
   const float x121 = Ts*F[52] + 1;
   const float x122 = Ts*F[50];
   const float x123 = x122*P[71] + P[67];
   const float x124 = Ts*F[53] + 1;
   const float x125 = Ts*F[51];
   const float x126 = x125*P[72] + P[69];
   const float x127 = Ts*F[54] + 1;
   Pnew[0] = x0*x1 + x0*P[3] + P[0] + Q[0];
   Pnew[1] = x1*x4 + x10*x9 + x11*x12 + x13*x14 + x15*x16 + x2*x3 + x5*x6 + x7*x8;
   Pnew[2] = x1*x18 + x10*x21 + x12*x22 + x14*x23 + x16*x17 + x19*x6 + x20*x8 + x24*x3;
   Pnew[3] = x1*x32 + x10*x29 + x12*x30 + x14*x31 + x16*x25 + x26*x3 + x27*x6 + x28*x8;
   Pnew[4] = x10*x34 + x12*x35 + x33*x8 + x36*x37 + x38*x39 + x40*x41 + x6;
   Pnew[5] = x10*x43 + x12*x44 + x37*x45 + x39*x46 + x41*x47 + x42*x6 + x8;
   Pnew[6] = x10 + x12*x50 + x37*x51 + x39*x52 + x41*x53 + x48*x6 + x49*x8;
   Pnew[7] = x10*x56 + x12 + x37*x57 + x39*x58 + x41*x59 + x54*x6 + x55*x8;
   Pnew[8] = x37;
   Pnew[9] = x39;
   Pnew[10] = x41;
   Pnew[11] = x14*x60;
   Pnew[12] = x11*x65 + x13*x61 + x15*x68 + x2*x66 + x4*x67 + x5*x62 + x63*x7 + x64*x9 + Q[1];
   Pnew[13] = x17*x68 + x18*x67 + x19*x62 + x20*x63 + x21*x64 + x22*x65 + x23*x61 + x24*x66;
   Pnew[14] = x25*x68 + x26*x66 + x27*x62 + x28*x63 + x29*x64 + x30*x65 + x31*x61 + x32*x67;
   Pnew[15] = x33*x63 + x34*x64 + x35*x65 + x36*x69 + x38*x70 + x40*x71 + x62;
   Pnew[16] = x42*x62 + x43*x64 + x44*x65 + x45*x69 + x46*x70 + x47*x71 + x63;
   Pnew[17] = x48*x62 + x49*x63 + x50*x65 + x51*x69 + x52*x70 + x53*x71 + x64;
   Pnew[18] = x54*x62 + x55*x63 + x56*x64 + x57*x69 + x58*x70 + x59*x71 + x65;
   Pnew[19] = x69;
   Pnew[20] = x70;
   Pnew[21] = x71;
   Pnew[22] = x60*x61;
   Pnew[23] = x17*x77 + x18*x78 + x19*x73 + x20*x74 + x21*x75 + x22*x76 + x23*x72 + x24*x79 + Q[2];
   Pnew[24] = x25*x77 + x26*x79 + x27*x73 + x28*x74 + x29*x75 + x30*x76 + x31*x72 + x32*x78;
   Pnew[25] = x33*x74 + x34*x75 + x35*x76 + x36*x80 + x38*x81 + x40*x82 + x73;
   Pnew[26] = x42*x73 + x43*x75 + x44*x76 + x45*x80 + x46*x81 + x47*x82 + x74;
   Pnew[27] = x48*x73 + x49*x74 + x50*x76 + x51*x80 + x52*x81 + x53*x82 + x75;
   Pnew[28] = x54*x73 + x55*x74 + x56*x75 + x57*x80 + x58*x81 + x59*x82 + x76;
   Pnew[29] = x80;
   Pnew[30] = x81;
   Pnew[31] = x82;
   Pnew[32] = x60*x72;
   Pnew[33] = x25*(x25*P[12] + x26*P[13] + x27*P[15] + x28*P[16] + x29*P[17] + x30*P[18] + x31*P[22] + x32*P[14]) + x26*(x25*P[13] + x26*P[23] + x27*P[25] + x28*P[26] + x29*P[27] + x30*P[28] + x31*P[32] + x32*P[24]) + x27*x84 + x28*x85 + x29*x86 + x30*x87 + x31*x83 + x32*(x25*P[14] + x26*P[24] + x27*P[34] + x28*P[35] + x29*P[36] + x30*P[37] + x31*P[41] + x32*P[33]) + Q[3];
   Pnew[34] = x33*x85 + x34*x86 + x35*x87 + x36*x88 + x38*x89 + x40*x90 + x84;
   Pnew[35] = x42*x84 + x43*x86 + x44*x87 + x45*x88 + x46*x89 + x47*x90 + x85;
   Pnew[36] = x48*x84 + x49*x85 + x50*x87 + x51*x88 + x52*x89 + x53*x90 + x86;
   Pnew[37] = x54*x84 + x55*x85 + x56*x86 + x57*x88 + x58*x89 + x59*x90 + x87;
   Pnew[38] = x88;
   Pnew[39] = x89;
   Pnew[40] = x90;
   Pnew[41] = x60*x83;
   Pnew[42] = x33*x95 + x34*x96 + x35*x97 + x36*x92 + x38*x93 + x40*x94 + x91 + Q[4];
   Pnew[43] = x42*x91 + x43*x96 + x44*x97 + x45*x92 + x46*x93 + x47*x94 + x95;
   Pnew[44] = x48*x91 + x49*x95 + x50*x97 + x51*x92 + x52*x93 + x53*x94 + x96;
   Pnew[45] = x54*x91 + x55*x95 + x56*x96 + x57*x92 + x58*x93 + x59*x94 + x97;
   Pnew[46] = x92 + x99*F[28];
   Pnew[47] = x100*F[29] + x93;
   Pnew[48] = x101*F[30] + x94;
   Pnew[49] = x102 + x103*x45 + x104*x46 + x105*x47 + x106*x42 + x107*x43 + x108*x44 + Q[5];
   Pnew[50] = x102*x49 + x103*x51 + x104*x52 + x105*x53 + x106*x48 + x107 + x108*x50;
   Pnew[51] = x102*x55 + x103*x57 + x104*x58 + x105*x59 + x106*x54 + x107*x56 + x108;
   Pnew[52] = x103 + x99*F[34];
   Pnew[53] = x100*F[35] + x104;
   Pnew[54] = x101*F[36] + x105;
   Pnew[55] = x109 + x110*x51 + x111*x52 + x112*x53 + x113*x48 + x114*x49 + x115*x50 + Q[6];
   Pnew[56] = x109*x56 + x110*x57 + x111*x58 + x112*x59 + x113*x54 + x114*x55 + x115;
   Pnew[57] = x110 + x99*F[40];
   Pnew[58] = x100*F[41] + x111;
   Pnew[59] = x101*F[42] + x112;
   Pnew[60] = x116*x57 + x117*x58 + x118*x59 + x54*(x54*P[42] + x55*P[43] + x56*P[44] + x57*P[46] + x58*P[47] + x59*P[48] + P[45]) + x54*P[45] + x55*(x54*P[43] + x55*P[49] + x56*P[50] + x57*P[52] + x58*P[53] + x59*P[54] + P[51]) + x55*P[51] + x56*(x54*P[44] + x55*P[50] + x56*P[55] + x57*P[57] + x58*P[58] + x59*P[59] + P[56]) + x56*P[56] + x57*P[61] + x58*P[62] + x59*P[63] + P[60] + Q[7];
   Pnew[61] = x116 + x99*F[46];
   Pnew[62] = x100*F[47] + x117;
   Pnew[63] = x101*F[48] + x118;
   Pnew[64] = x119*x120 + x119*P[65] + P[64] + Q[8];
   Pnew[65] = x120*x121;
   Pnew[66] = x122*x123 + x122*P[67] + P[66] + Q[9];
   Pnew[67] = x123*x124;
   Pnew[68] = x125*x126 + x125*P[69] + P[68] + Q[10];
   Pnew[69] = x126*x127;
   Pnew[70] = (x121*x121)*P[70] + Q[11];
   Pnew[71] = (x124*x124)*P[71] + Q[12];
   Pnew[72] = (x127*x127)*P[72] + Q[13];
   Pnew[73] = (x60*x60)*P[73] + Q[14];
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
   Pnew[49] = -x5*(P[5]*P[5]) + P[49];
   Pnew[50] = -x10*P[6] + P[50];
   Pnew[51] = -x10*P[7] + P[51];
   Pnew[52] = -x10*P[8] + P[52];
   Pnew[53] = -x10*P[9] + P[53];
   Pnew[54] = -x10*P[10] + P[54];
   Pnew[55] = -x5*(P[6]*P[6]) + P[55];
   Pnew[56] = -x11*P[7] + P[56];
   Pnew[57] = -x11*P[8] + P[57];
   Pnew[58] = -x11*P[9] + P[58];
   Pnew[59] = -x11*P[10] + P[59];
   Pnew[60] = -x5*(P[7]*P[7]) + P[60];
   Pnew[61] = -x12*P[8] + P[61];
   Pnew[62] = -x12*P[9] + P[62];
   Pnew[63] = -x12*P[10] + P[63];
   Pnew[64] = -x5*(P[8]*P[8]) + P[64];
   Pnew[65] = P[65];
   Pnew[66] = -x5*(P[9]*P[9]) + P[66];
   Pnew[67] = P[67];
   Pnew[68] = -x5*(P[10]*P[10]) + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = P[72];
   Pnew[73] = -x5*(P[11]*P[11]) + P[73];
}

void mag_correction(const float *restrict x, const float *restrict P, const float *restrict mag,
        const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = -mag[1] - 2*x[4]*x[7] + 2*x[5]*x[6];
   const float x1 = H[19]*P[4] + H[20]*P[5] + H[21]*P[6] + H[22]*P[7];
   const float x2 = H[19]*P[42] + H[20]*P[43] + H[21]*P[44] + H[22]*P[45];
   const float x3 = x2*H[23];
   const float x4 = H[19]*P[43] + H[20]*P[49] + H[21]*P[50] + H[22]*P[51];
   const float x5 = x4*H[24];
   const float x6 = H[19]*P[44] + H[20]*P[50] + H[21]*P[55] + H[22]*P[56];
   const float x7 = x6*H[25];
   const float x8 = H[19]*P[45] + H[20]*P[51] + H[21]*P[56] + H[22]*P[60];
   const float x9 = x8*H[26];
   const float x10 = H[23]*P[42] + H[24]*P[43] + H[25]*P[44] + H[26]*P[45];
   const float x11 = x10*H[19];
   const float x12 = H[23]*P[43] + H[24]*P[49] + H[25]*P[50] + H[26]*P[51];
   const float x13 = x12*H[20];
   const float x14 = H[23]*P[44] + H[24]*P[50] + H[25]*P[55] + H[26]*P[56];
   const float x15 = x14*H[21];
   const float x16 = H[23]*P[45] + H[24]*P[51] + H[25]*P[56] + H[26]*P[60];
   const float x17 = x16*H[22];
   const float x18 = x2*H[19] + x4*H[20] + x6*H[21] + x8*H[22] + R[7];
   const float x19 = x10*H[23] + x12*H[24] + x14*H[25] + x16*H[26] + R[8];
   const float x20 = 1.0F/(x18*x19 - (x11 + x13 + x15 + x17)*(x3 + x5 + x7 + x9));
   const float x21 = x20*(-x3 - x5 - x7 - x9);
   const float x22 = x1*x21;
   const float x23 = H[23]*P[4] + H[24]*P[5] + H[25]*P[6] + H[26]*P[7];
   const float x24 = x18*x20;
   const float x25 = x23*x24;
   const float x26 = -mag[0] + (x[4]*x[4]) + (x[5]*x[5]) - (x[6]*x[6]) - (x[7]*x[7]);
   const float x27 = x20*(-x11 - x13 - x15 - x17);
   const float x28 = x23*x27;
   const float x29 = x19*x20;
   const float x30 = x1*x29;
   const float x31 = H[19]*P[15] + H[20]*P[16] + H[21]*P[17] + H[22]*P[18];
   const float x32 = x21*x31;
   const float x33 = H[23]*P[15] + H[24]*P[16] + H[25]*P[17] + H[26]*P[18];
   const float x34 = x24*x33;
   const float x35 = x27*x33;
   const float x36 = x29*x31;
   const float x37 = H[19]*P[25] + H[20]*P[26] + H[21]*P[27] + H[22]*P[28];
   const float x38 = x21*x37;
   const float x39 = H[23]*P[25] + H[24]*P[26] + H[25]*P[27] + H[26]*P[28];
   const float x40 = x24*x39;
   const float x41 = x27*x39;
   const float x42 = x29*x37;
   const float x43 = H[19]*P[34] + H[20]*P[35] + H[21]*P[36] + H[22]*P[37];
   const float x44 = x21*x43;
   const float x45 = H[23]*P[34] + H[24]*P[35] + H[25]*P[36] + H[26]*P[37];
   const float x46 = x24*x45;
   const float x47 = x27*x45;
   const float x48 = x29*x43;
   const float x49 = x2*x21;
   const float x50 = x10*x24;
   const float x51 = x10*x27;
   const float x52 = x2*x29;
   const float x53 = x21*x4;
   const float x54 = x12*x24;
   const float x55 = x12*x27;
   const float x56 = x29*x4;
   const float x57 = x21*x6;
   const float x58 = x14*x24;
   const float x59 = x14*x27;
   const float x60 = x29*x6;
   const float x61 = x21*x8;
   const float x62 = x16*x24;
   const float x63 = x16*x27;
   const float x64 = x29*x8;
   const float x65 = H[19]*P[46] + H[20]*P[52] + H[21]*P[57] + H[22]*P[61];
   const float x66 = x21*x65;
   const float x67 = H[23]*P[46] + H[24]*P[52] + H[25]*P[57] + H[26]*P[61];
   const float x68 = x24*x67;
   const float x69 = x27*x67;
   const float x70 = x29*x65;
   const float x71 = H[19]*P[47] + H[20]*P[53] + H[21]*P[58] + H[22]*P[62];
   const float x72 = x21*x71;
   const float x73 = H[23]*P[47] + H[24]*P[53] + H[25]*P[58] + H[26]*P[62];
   const float x74 = x24*x73;
   const float x75 = x27*x73;
   const float x76 = x29*x71;
   const float x77 = H[19]*P[48] + H[20]*P[54] + H[21]*P[59] + H[22]*P[63];
   const float x78 = x21*x77;
   const float x79 = H[23]*P[48] + H[24]*P[54] + H[25]*P[59] + H[26]*P[63];
   const float x80 = x24*x79;
   const float x81 = x27*x79;
   const float x82 = x29*x77;
   const float x83 = x22 + x25;
   const float x84 = x28 + x30;
   const float x85 = -x83*H[23] - x84*H[19];
   const float x86 = -x83*H[24] - x84*H[20];
   const float x87 = -x83*H[25] - x84*H[21];
   const float x88 = -x83*H[26] - x84*H[22];
   const float x89 = x32 + x34;
   const float x90 = x35 + x36;
   const float x91 = -x89*H[23] - x90*H[19];
   const float x92 = -x89*H[24] - x90*H[20];
   const float x93 = -x89*H[25] - x90*H[21];
   const float x94 = -x89*H[26] - x90*H[22];
   const float x95 = x38 + x40;
   const float x96 = x41 + x42;
   const float x97 = -x95*H[23] - x96*H[19];
   const float x98 = -x95*H[24] - x96*H[20];
   const float x99 = -x95*H[25] - x96*H[21];
   const float x100 = -x95*H[26] - x96*H[22];
   const float x101 = x44 + x46;
   const float x102 = x47 + x48;
   const float x103 = -x101*H[23] - x102*H[19];
   const float x104 = -x101*H[24] - x102*H[20];
   const float x105 = -x101*H[25] - x102*H[21];
   const float x106 = -x101*H[26] - x102*H[22];
   const float x107 = x49 + x50;
   const float x108 = x51 + x52;
   const float x109 = -x107*H[24] - x108*H[20];
   const float x110 = -x107*H[25] - x108*H[21];
   const float x111 = -x107*H[26] - x108*H[22];
   const float x112 = -x107*H[23] - x108*H[19] + 1;
   const float x113 = x53 + x54;
   const float x114 = x55 + x56;
   const float x115 = -x113*H[23] - x114*H[19];
   const float x116 = -x113*H[25] - x114*H[21];
   const float x117 = -x113*H[26] - x114*H[22];
   const float x118 = -x113*H[24] - x114*H[20] + 1;
   const float x119 = x57 + x58;
   const float x120 = x59 + x60;
   const float x121 = -x119*H[23] - x120*H[19];
   const float x122 = -x119*H[24] - x120*H[20];
   const float x123 = -x119*H[26] - x120*H[22];
   const float x124 = -x119*H[25] - x120*H[21] + 1;
   const float x125 = x61 + x62;
   const float x126 = x63 + x64;
   const float x127 = -x125*H[23] - x126*H[19];
   const float x128 = -x125*H[24] - x126*H[20];
   const float x129 = -x125*H[25] - x126*H[21];
   const float x130 = -x125*H[26] - x126*H[22] + 1;
   const float x131 = x66 + x68;
   const float x132 = x69 + x70;
   const float x133 = x72 + x74;
   const float x134 = x75 + x76;
   const float x135 = x78 + x80;
   const float x136 = x81 + x82;
   xnew[0] = x0*(-x22 - x25) + x26*(-x28 - x30) + x[0];
   xnew[1] = x0*(-x32 - x34) + x26*(-x35 - x36) + x[1];
   xnew[2] = x0*(-x38 - x40) + x26*(-x41 - x42) + x[2];
   xnew[3] = x0*(-x44 - x46) + x26*(-x47 - x48) + x[3];
   xnew[4] = x0*(-x49 - x50) + x26*(-x51 - x52) + x[4];
   xnew[5] = x0*(-x53 - x54) + x26*(-x55 - x56) + x[5];
   xnew[6] = x0*(-x57 - x58) + x26*(-x59 - x60) + x[6];
   xnew[7] = x0*(-x61 - x62) + x26*(-x63 - x64) + x[7];
   xnew[8] = x0*(-x66 - x68) + x26*(-x69 - x70) + x[8];
   xnew[9] = x0*(-x72 - x74) + x26*(-x75 - x76) + x[9];
   xnew[10] = x0*(-x78 - x80) + x26*(-x81 - x82) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x[14];
   Pnew[0] = x85*P[4] + x86*P[5] + x87*P[6] + x88*P[7] + P[0];
   Pnew[1] = x85*P[15] + x86*P[16] + x87*P[17] + x88*P[18] + P[1];
   Pnew[2] = x85*P[25] + x86*P[26] + x87*P[27] + x88*P[28] + P[2];
   Pnew[3] = x85*P[34] + x86*P[35] + x87*P[36] + x88*P[37] + P[3];
   Pnew[4] = x85*P[42] + x86*P[43] + x87*P[44] + x88*P[45] + P[4];
   Pnew[5] = x85*P[43] + x86*P[49] + x87*P[50] + x88*P[51] + P[5];
   Pnew[6] = x85*P[44] + x86*P[50] + x87*P[55] + x88*P[56] + P[6];
   Pnew[7] = x85*P[45] + x86*P[51] + x87*P[56] + x88*P[60] + P[7];
   Pnew[8] = x85*P[46] + x86*P[52] + x87*P[57] + x88*P[61] + P[8];
   Pnew[9] = x85*P[47] + x86*P[53] + x87*P[58] + x88*P[62] + P[9];
   Pnew[10] = x85*P[48] + x86*P[54] + x87*P[59] + x88*P[63] + P[10];
   Pnew[11] = P[11];
   Pnew[12] = x91*P[15] + x92*P[16] + x93*P[17] + x94*P[18] + P[12];
   Pnew[13] = x91*P[25] + x92*P[26] + x93*P[27] + x94*P[28] + P[13];
   Pnew[14] = x91*P[34] + x92*P[35] + x93*P[36] + x94*P[37] + P[14];
   Pnew[15] = x91*P[42] + x92*P[43] + x93*P[44] + x94*P[45] + P[15];
   Pnew[16] = x91*P[43] + x92*P[49] + x93*P[50] + x94*P[51] + P[16];
   Pnew[17] = x91*P[44] + x92*P[50] + x93*P[55] + x94*P[56] + P[17];
   Pnew[18] = x91*P[45] + x92*P[51] + x93*P[56] + x94*P[60] + P[18];
   Pnew[19] = x91*P[46] + x92*P[52] + x93*P[57] + x94*P[61] + P[19];
   Pnew[20] = x91*P[47] + x92*P[53] + x93*P[58] + x94*P[62] + P[20];
   Pnew[21] = x91*P[48] + x92*P[54] + x93*P[59] + x94*P[63] + P[21];
   Pnew[22] = P[22];
   Pnew[23] = x100*P[28] + x97*P[25] + x98*P[26] + x99*P[27] + P[23];
   Pnew[24] = x100*P[37] + x97*P[34] + x98*P[35] + x99*P[36] + P[24];
   Pnew[25] = x100*P[45] + x97*P[42] + x98*P[43] + x99*P[44] + P[25];
   Pnew[26] = x100*P[51] + x97*P[43] + x98*P[49] + x99*P[50] + P[26];
   Pnew[27] = x100*P[56] + x97*P[44] + x98*P[50] + x99*P[55] + P[27];
   Pnew[28] = x100*P[60] + x97*P[45] + x98*P[51] + x99*P[56] + P[28];
   Pnew[29] = x100*P[61] + x97*P[46] + x98*P[52] + x99*P[57] + P[29];
   Pnew[30] = x100*P[62] + x97*P[47] + x98*P[53] + x99*P[58] + P[30];
   Pnew[31] = x100*P[63] + x97*P[48] + x98*P[54] + x99*P[59] + P[31];
   Pnew[32] = P[32];
   Pnew[33] = x103*P[34] + x104*P[35] + x105*P[36] + x106*P[37] + P[33];
   Pnew[34] = x103*P[42] + x104*P[43] + x105*P[44] + x106*P[45] + P[34];
   Pnew[35] = x103*P[43] + x104*P[49] + x105*P[50] + x106*P[51] + P[35];
   Pnew[36] = x103*P[44] + x104*P[50] + x105*P[55] + x106*P[56] + P[36];
   Pnew[37] = x103*P[45] + x104*P[51] + x105*P[56] + x106*P[60] + P[37];
   Pnew[38] = x103*P[46] + x104*P[52] + x105*P[57] + x106*P[61] + P[38];
   Pnew[39] = x103*P[47] + x104*P[53] + x105*P[58] + x106*P[62] + P[39];
   Pnew[40] = x103*P[48] + x104*P[54] + x105*P[59] + x106*P[63] + P[40];
   Pnew[41] = P[41];
   Pnew[42] = x109*P[43] + x110*P[44] + x111*P[45] + x112*P[42];
   Pnew[43] = x109*P[49] + x110*P[50] + x111*P[51] + x112*P[43];
   Pnew[44] = x109*P[50] + x110*P[55] + x111*P[56] + x112*P[44];
   Pnew[45] = x109*P[51] + x110*P[56] + x111*P[60] + x112*P[45];
   Pnew[46] = x109*P[52] + x110*P[57] + x111*P[61] + x112*P[46];
   Pnew[47] = x109*P[53] + x110*P[58] + x111*P[62] + x112*P[47];
   Pnew[48] = x109*P[54] + x110*P[59] + x111*P[63] + x112*P[48];
   Pnew[49] = x115*P[43] + x116*P[50] + x117*P[51] + x118*P[49];
   Pnew[50] = x115*P[44] + x116*P[55] + x117*P[56] + x118*P[50];
   Pnew[51] = x115*P[45] + x116*P[56] + x117*P[60] + x118*P[51];
   Pnew[52] = x115*P[46] + x116*P[57] + x117*P[61] + x118*P[52];
   Pnew[53] = x115*P[47] + x116*P[58] + x117*P[62] + x118*P[53];
   Pnew[54] = x115*P[48] + x116*P[59] + x117*P[63] + x118*P[54];
   Pnew[55] = x121*P[44] + x122*P[50] + x123*P[56] + x124*P[55];
   Pnew[56] = x121*P[45] + x122*P[51] + x123*P[60] + x124*P[56];
   Pnew[57] = x121*P[46] + x122*P[52] + x123*P[61] + x124*P[57];
   Pnew[58] = x121*P[47] + x122*P[53] + x123*P[62] + x124*P[58];
   Pnew[59] = x121*P[48] + x122*P[54] + x123*P[63] + x124*P[59];
   Pnew[60] = x127*P[45] + x128*P[51] + x129*P[56] + x130*P[60];
   Pnew[61] = x127*P[46] + x128*P[52] + x129*P[57] + x130*P[61];
   Pnew[62] = x127*P[47] + x128*P[53] + x129*P[58] + x130*P[62];
   Pnew[63] = x127*P[48] + x128*P[54] + x129*P[59] + x130*P[63];
   Pnew[64] = (-x131*H[23] - x132*H[19])*P[46] + (-x131*H[24] - x132*H[20])*P[52] + (-x131*H[25] - x132*H[21])*P[57] + (-x131*H[26] - x132*H[22])*P[61] + P[64];
   Pnew[65] = P[65];
   Pnew[66] = (-x133*H[23] - x134*H[19])*P[47] + (-x133*H[24] - x134*H[20])*P[53] + (-x133*H[25] - x134*H[21])*P[58] + (-x133*H[26] - x134*H[22])*P[62] + P[66];
   Pnew[67] = P[67];
   Pnew[68] = (-x135*H[23] - x136*H[19])*P[48] + (-x135*H[24] - x136*H[20])*P[54] + (-x135*H[25] - x136*H[21])*P[59] + (-x135*H[26] - x136*H[22])*P[63] + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = P[72];
   Pnew[73] = P[73];
}

void gyro_correction(const float * restrict x, const float * restrict P, const float * restrict gyro,
	const float *restrict H, const float *restrict R, float *restrict xnew, float *restrict Pnew)
{
   const float x0 = (H[16]*H[16]);
   const float x1 = x0*P[64];
   const float x2 = 1.0F/(x1 + R[4]);
   const float x3 = x2*(-gyro[0] + x[8])*H[16];
   const float x4 = (H[17]*H[17]);
   const float x5 = x4*P[66];
   const float x6 = 1.0F/(x5 + R[5]);
   const float x7 = x6*(-gyro[1] + x[9])*H[17];
   const float x8 = (H[18]*H[18]);
   const float x9 = x8*P[68];
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
   const float x33 = x0*x2*P[52];
   const float x34 = x4*x6*P[53];
   const float x35 = x10*x8*P[54];
   const float x36 = -x18 + 1;
   const float x37 = -x19 + 1;
   const float x38 = -x20 + 1;
   xnew[0] = -x11*P[10] - x3*P[8] - x7*P[9] + x[0];
   xnew[1] = -x11*P[21] - x3*P[19] - x7*P[20] + x[1];
   xnew[2] = -x11*P[31] - x3*P[29] - x7*P[30] + x[2];
   xnew[3] = -x11*P[40] - x3*P[38] - x7*P[39] + x[3];
   xnew[4] = -x11*P[48] - x3*P[46] - x7*P[47] + x[4];
   xnew[5] = -x11*P[54] - x3*P[52] - x7*P[53] + x[5];
   xnew[6] = -x11*P[59] - x3*P[57] - x7*P[58] + x[6];
   xnew[7] = -x11*P[63] - x3*P[61] - x7*P[62] + x[7];
   xnew[8] = -x3*P[64] + x[8];
   xnew[9] = -x7*P[66] + x[9];
   xnew[10] = -x11*P[68] + x[10];
   xnew[11] = -x3*P[65] + x[11];
   xnew[12] = -x7*P[67] + x[12];
   xnew[13] = -x11*P[69] + x[13];
   xnew[14] = x[14];
   Pnew[0] = -x12*(P[8]*P[8]) - x13*(P[9]*P[9]) - x14*(P[10]*P[10]) + P[0];
   Pnew[1] = -x15*P[19] - x16*P[20] - x17*P[21] + P[1];
   Pnew[2] = -x15*P[29] - x16*P[30] - x17*P[31] + P[2];
   Pnew[3] = -x15*P[38] - x16*P[39] - x17*P[40] + P[3];
   Pnew[4] = -x15*P[46] - x16*P[47] - x17*P[48] + P[4];
   Pnew[5] = -x15*P[52] - x16*P[53] - x17*P[54] + P[5];
   Pnew[6] = -x15*P[57] - x16*P[58] - x17*P[59] + P[6];
   Pnew[7] = -x15*P[61] - x16*P[62] - x17*P[63] + P[7];
   Pnew[8] = -x18*P[8] + P[8];
   Pnew[9] = -x19*P[9] + P[9];
   Pnew[10] = -x20*P[10] + P[10];
   Pnew[11] = P[11];
   Pnew[12] = -x12*(P[19]*P[19]) - x13*(P[20]*P[20]) - x14*(P[21]*P[21]) + P[12];
   Pnew[13] = -x21*P[29] - x22*P[30] - x23*P[31] + P[13];
   Pnew[14] = -x21*P[38] - x22*P[39] - x23*P[40] + P[14];
   Pnew[15] = -x21*P[46] - x22*P[47] - x23*P[48] + P[15];
   Pnew[16] = -x21*P[52] - x22*P[53] - x23*P[54] + P[16];
   Pnew[17] = -x21*P[57] - x22*P[58] - x23*P[59] + P[17];
   Pnew[18] = -x21*P[61] - x22*P[62] - x23*P[63] + P[18];
   Pnew[19] = -x18*P[19] + P[19];
   Pnew[20] = -x19*P[20] + P[20];
   Pnew[21] = -x20*P[21] + P[21];
   Pnew[22] = P[22];
   Pnew[23] = -x12*(P[29]*P[29]) - x13*(P[30]*P[30]) - x14*(P[31]*P[31]) + P[23];
   Pnew[24] = -x24*P[38] - x25*P[39] - x26*P[40] + P[24];
   Pnew[25] = -x24*P[46] - x25*P[47] - x26*P[48] + P[25];
   Pnew[26] = -x24*P[52] - x25*P[53] - x26*P[54] + P[26];
   Pnew[27] = -x24*P[57] - x25*P[58] - x26*P[59] + P[27];
   Pnew[28] = -x24*P[61] - x25*P[62] - x26*P[63] + P[28];
   Pnew[29] = -x18*P[29] + P[29];
   Pnew[30] = -x19*P[30] + P[30];
   Pnew[31] = -x20*P[31] + P[31];
   Pnew[32] = P[32];
   Pnew[33] = -x12*(P[38]*P[38]) - x13*(P[39]*P[39]) - x14*(P[40]*P[40]) + P[33];
   Pnew[34] = -x27*P[46] - x28*P[47] - x29*P[48] + P[34];
   Pnew[35] = -x27*P[52] - x28*P[53] - x29*P[54] + P[35];
   Pnew[36] = -x27*P[57] - x28*P[58] - x29*P[59] + P[36];
   Pnew[37] = -x27*P[61] - x28*P[62] - x29*P[63] + P[37];
   Pnew[38] = -x18*P[38] + P[38];
   Pnew[39] = -x19*P[39] + P[39];
   Pnew[40] = -x20*P[40] + P[40];
   Pnew[41] = P[41];
   Pnew[42] = -x12*(P[46]*P[46]) - x13*(P[47]*P[47]) - x14*(P[48]*P[48]) + P[42];
   Pnew[43] = -x30*P[52] - x31*P[53] - x32*P[54] + P[43];
   Pnew[44] = -x30*P[57] - x31*P[58] - x32*P[59] + P[44];
   Pnew[45] = -x30*P[61] - x31*P[62] - x32*P[63] + P[45];
   Pnew[46] = -x18*P[46] + P[46];
   Pnew[47] = -x19*P[47] + P[47];
   Pnew[48] = -x20*P[48] + P[48];
   Pnew[49] = -x12*(P[52]*P[52]) - x13*(P[53]*P[53]) - x14*(P[54]*P[54]) + P[49];
   Pnew[50] = -x33*P[57] - x34*P[58] - x35*P[59] + P[50];
   Pnew[51] = -x33*P[61] - x34*P[62] - x35*P[63] + P[51];
   Pnew[52] = -x18*P[52] + P[52];
   Pnew[53] = -x19*P[53] + P[53];
   Pnew[54] = -x20*P[54] + P[54];
   Pnew[55] = -x12*(P[57]*P[57]) - x13*(P[58]*P[58]) - x14*(P[59]*P[59]) + P[55];
   Pnew[56] = -x12*P[57]*P[61] - x13*P[58]*P[62] - x14*P[59]*P[63] + P[56];
   Pnew[57] = -x18*P[57] + P[57];
   Pnew[58] = -x19*P[58] + P[58];
   Pnew[59] = -x20*P[59] + P[59];
   Pnew[60] = -x12*(P[61]*P[61]) - x13*(P[62]*P[62]) - x14*(P[63]*P[63]) + P[60];
   Pnew[61] = -x18*P[61] + P[61];
   Pnew[62] = -x19*P[62] + P[62];
   Pnew[63] = -x20*P[63] + P[63];
   Pnew[64] = x36*P[64];
   Pnew[65] = x36*P[65];
   Pnew[66] = x37*P[66];
   Pnew[67] = x37*P[67];
   Pnew[68] = x38*P[68];
   Pnew[69] = x38*P[69];
   Pnew[70] = -x12*(P[65]*P[65]) + P[70];
   Pnew[71] = -x13*(P[67]*P[67]) + P[71];
   Pnew[72] = -x14*(P[69]*P[69]) + P[72];
   Pnew[73] = P[73];
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
   const float x6 = -(x0*(x[4]*x[5] + x[6]*x[7]) - 2*(x1 - x2)*x[1] + (x3 + x4 - x5)*x[2])*param[3] - accel[1];
   const float x7 = H[1]*P[1] + H[2]*P[2] + H[3]*P[3] + H[4]*P[4] + H[5]*P[5] + H[6]*P[6] + H[7]*P[7];
   const float x8 = H[1]*P[22];
   const float x9 = H[2]*P[32];
   const float x10 = H[3]*P[41];
   const float x11 = x10 + x8 + x9;
   const float x12 = x11*H[15];
   const float x13 = H[8]*P[22];
   const float x14 = H[9]*P[32];
   const float x15 = H[10]*P[41];
   const float x16 = x13*H[15] + x14*H[15] + x15*H[15];
   const float x17 = x12*x16;
   const float x18 = (H[15]*H[15])*P[73] + R[3];
   const float x19 = H[1]*P[12] + H[2]*P[13] + H[3]*P[14] + H[4]*P[15] + H[5]*P[16] + H[6]*P[17] + H[7]*P[18];
   const float x20 = H[1]*P[13] + H[2]*P[23] + H[3]*P[24] + H[4]*P[25] + H[5]*P[26] + H[6]*P[27] + H[7]*P[28];
   const float x21 = H[1]*P[14] + H[2]*P[24] + H[3]*P[33] + H[4]*P[34] + H[5]*P[35] + H[6]*P[36] + H[7]*P[37];
   const float x22 = H[1]*P[15] + H[2]*P[25] + H[3]*P[34] + H[4]*P[42] + H[5]*P[43] + H[6]*P[44] + H[7]*P[45];
   const float x23 = H[1]*P[16] + H[2]*P[26] + H[3]*P[35] + H[4]*P[43] + H[5]*P[49] + H[6]*P[50] + H[7]*P[51];
   const float x24 = H[1]*P[17] + H[2]*P[27] + H[3]*P[36] + H[4]*P[44] + H[5]*P[50] + H[6]*P[55] + H[7]*P[56];
   const float x25 = H[1]*P[18] + H[2]*P[28] + H[3]*P[37] + H[4]*P[45] + H[5]*P[51] + H[6]*P[56] + H[7]*P[60];
   const float x26 = x19*H[8] + x20*H[9] + x21*H[10] + x22*H[11] + x23*H[12] + x24*H[13] + x25*H[14];
   const float x27 = x18*x26;
   const float x28 = x17 - x27;
   const float x29 = H[8]*P[12] + H[9]*P[13] + H[10]*P[14] + H[11]*P[15] + H[12]*P[16] + H[13]*P[17] + H[14]*P[18];
   const float x30 = H[8]*P[13] + H[9]*P[23] + H[10]*P[24] + H[11]*P[25] + H[12]*P[26] + H[13]*P[27] + H[14]*P[28];
   const float x31 = H[8]*P[14] + H[9]*P[24] + H[10]*P[33] + H[11]*P[34] + H[12]*P[35] + H[13]*P[36] + H[14]*P[37];
   const float x32 = H[8]*P[15] + H[9]*P[25] + H[10]*P[34] + H[11]*P[42] + H[12]*P[43] + H[13]*P[44] + H[14]*P[45];
   const float x33 = H[8]*P[16] + H[9]*P[26] + H[10]*P[35] + H[11]*P[43] + H[12]*P[49] + H[13]*P[50] + H[14]*P[51];
   const float x34 = H[8]*P[17] + H[9]*P[27] + H[10]*P[36] + H[11]*P[44] + H[12]*P[50] + H[13]*P[55] + H[14]*P[56];
   const float x35 = H[8]*P[18] + H[9]*P[28] + H[10]*P[37] + H[11]*P[45] + H[12]*P[51] + H[13]*P[56] + H[14]*P[60];
   const float x36 = x29*H[1] + x30*H[2] + x31*H[3] + x32*H[4] + x33*H[5] + x34*H[6] + x35*H[7];
   const float x37 = x13 + x14 + x15;
   const float x38 = x37*H[15];
   const float x39 = x10*H[15] + x8*H[15] + x9*H[15];
   const float x40 = x38*x39;
   const float x41 = x12*x39;
   const float x42 = x29*H[8] + x30*H[9] + x31*H[10] + x32*H[11] + x33*H[12] + x34*H[13] + x35*H[14] + R[2];
   const float x43 = x16*x38;
   const float x44 = x19*H[1] + x20*H[2] + x21*H[3] + x22*H[4] + x23*H[5] + x24*H[6] + x25*H[7] + R[1];
   const float x45 = x18*x44;
   const float x46 = 1.0F/(x17*x36 + x26*x40 - x27*x36 - x41*x42 + x42*x45 - x43*x44);
   const float x47 = x28*x46;
   const float x48 = x47*x7;
   const float x49 = H[8]*P[1] + H[9]*P[2] + H[10]*P[3] + H[11]*P[4] + H[12]*P[5] + H[13]*P[6] + H[14]*P[7];
   const float x50 = -x41 + x45;
   const float x51 = x46*x50;
   const float x52 = x49*x51;
   const float x53 = x46*(-x16*x44 + x26*x39)*H[15];
   const float x54 = x53*P[11];
   const float x55 = -(-x0*(x[4]*x[6] - x[5]*x[7]) + 2*(x1 + x2)*x[2] + (x3 - x4 + x5)*x[1])*param[3] - accel[0];
   const float x56 = -x18*x36 + x40;
   const float x57 = x46*x56;
   const float x58 = x49*x57;
   const float x59 = x18*x42 - x43;
   const float x60 = x46*x59;
   const float x61 = x60*x7;
   const float x62 = x46*(x16*x36 - x39*x42)*H[15];
   const float x63 = x62*P[11];
   const float x64 = -accel[2] - param[0]*x[14];
   const float x65 = x12*x36 - x38*x44;
   const float x66 = x46*x65;
   const float x67 = x49*x66;
   const float x68 = -x12*x42 + x26*x38;
   const float x69 = x46*x68;
   const float x70 = x69*x7;
   const float x71 = x46*(-x26*x36 + x42*x44)*H[15];
   const float x72 = x71*P[11];
   const float x73 = x19*x47;
   const float x74 = x29*x51;
   const float x75 = x53*P[22];
   const float x76 = x29*x57;
   const float x77 = x19*x60;
   const float x78 = x62*P[22];
   const float x79 = x29*x66;
   const float x80 = x19*x69;
   const float x81 = x71*P[22];
   const float x82 = x20*x47;
   const float x83 = x30*x51;
   const float x84 = x53*P[32];
   const float x85 = x30*x57;
   const float x86 = x20*x60;
   const float x87 = x62*P[32];
   const float x88 = x30*x66;
   const float x89 = x20*x69;
   const float x90 = x71*P[32];
   const float x91 = x21*x47;
   const float x92 = x31*x51;
   const float x93 = x53*P[41];
   const float x94 = x31*x57;
   const float x95 = x21*x60;
   const float x96 = x62*P[41];
   const float x97 = x31*x66;
   const float x98 = x21*x69;
   const float x99 = x71*P[41];
   const float x100 = x22*x47;
   const float x101 = x32*x51;
   const float x102 = x32*x57;
   const float x103 = x22*x60;
   const float x104 = x23*x47;
   const float x105 = x33*x51;
   const float x106 = x33*x57;
   const float x107 = x23*x60;
   const float x108 = x24*x47;
   const float x109 = x34*x51;
   const float x110 = x34*x57;
   const float x111 = x24*x60;
   const float x112 = x25*x47;
   const float x113 = x35*x51;
   const float x114 = x35*x57;
   const float x115 = x25*x60;
   const float x116 = H[1]*P[19] + H[2]*P[29] + H[3]*P[38] + H[4]*P[46] + H[5]*P[52] + H[6]*P[57] + H[7]*P[61];
   const float x117 = x116*x47;
   const float x118 = H[8]*P[19] + H[9]*P[29] + H[10]*P[38] + H[11]*P[46] + H[12]*P[52] + H[13]*P[57] + H[14]*P[61];
   const float x119 = x118*x51;
   const float x120 = x118*x57;
   const float x121 = x116*x60;
   const float x122 = H[1]*P[20] + H[2]*P[30] + H[3]*P[39] + H[4]*P[47] + H[5]*P[53] + H[6]*P[58] + H[7]*P[62];
   const float x123 = x122*x47;
   const float x124 = H[8]*P[20] + H[9]*P[30] + H[10]*P[39] + H[11]*P[47] + H[12]*P[53] + H[13]*P[58] + H[14]*P[62];
   const float x125 = x124*x51;
   const float x126 = x124*x57;
   const float x127 = x122*x60;
   const float x128 = H[1]*P[21] + H[2]*P[31] + H[3]*P[40] + H[4]*P[48] + H[5]*P[54] + H[6]*P[59] + H[7]*P[63];
   const float x129 = x128*x47;
   const float x130 = H[8]*P[21] + H[9]*P[31] + H[10]*P[40] + H[11]*P[48] + H[12]*P[54] + H[13]*P[59] + H[14]*P[63];
   const float x131 = x130*x51;
   const float x132 = x130*x57;
   const float x133 = x128*x60;
   const float x134 = x11*x46;
   const float x135 = x134*x28;
   const float x136 = x37*x46;
   const float x137 = x136*x50;
   const float x138 = x53*P[73];
   const float x139 = x136*x56;
   const float x140 = x134*x59;
   const float x141 = x62*P[73];
   const float x142 = x136*x65;
   const float x143 = x134*x68;
   const float x144 = x71*P[73];
   const float x145 = (x67 + x70 + x72)*H[15];
   const float x146 = x48 + x52 + x54;
   const float x147 = x58 + x61 + x63;
   const float x148 = -x146*H[8] - x147*H[1];
   const float x149 = -x146*H[9] - x147*H[2];
   const float x150 = -x146*H[10] - x147*H[3];
   const float x151 = -x146*H[11] - x147*H[4];
   const float x152 = -x146*H[12] - x147*H[5];
   const float x153 = -x146*H[13] - x147*H[6];
   const float x154 = -x146*H[14] - x147*H[7];
   const float x155 = (x79 + x80 + x81)*H[15];
   const float x156 = x73 + x74 + x75;
   const float x157 = x76 + x77 + x78;
   const float x158 = -x156*H[9] - x157*H[2];
   const float x159 = -x156*H[10] - x157*H[3];
   const float x160 = -x156*H[11] - x157*H[4];
   const float x161 = -x156*H[12] - x157*H[5];
   const float x162 = -x156*H[13] - x157*H[6];
   const float x163 = -x156*H[14] - x157*H[7];
   const float x164 = -x156*H[8] - x157*H[1] + 1;
   const float x165 = (x88 + x89 + x90)*H[15];
   const float x166 = x82 + x83 + x84;
   const float x167 = x85 + x86 + x87;
   const float x168 = -x166*H[8] - x167*H[1];
   const float x169 = -x166*H[10] - x167*H[3];
   const float x170 = -x166*H[11] - x167*H[4];
   const float x171 = -x166*H[12] - x167*H[5];
   const float x172 = -x166*H[13] - x167*H[6];
   const float x173 = -x166*H[14] - x167*H[7];
   const float x174 = -x166*H[9] - x167*H[2] + 1;
   const float x175 = (x97 + x98 + x99)*H[15];
   const float x176 = x91 + x92 + x93;
   const float x177 = x94 + x95 + x96;
   const float x178 = -x176*H[8] - x177*H[1];
   const float x179 = -x176*H[9] - x177*H[2];
   const float x180 = -x176*H[11] - x177*H[4];
   const float x181 = -x176*H[12] - x177*H[5];
   const float x182 = -x176*H[13] - x177*H[6];
   const float x183 = -x176*H[14] - x177*H[7];
   const float x184 = -x176*H[10] - x177*H[3] + 1;
   const float x185 = x100 + x101;
   const float x186 = x102 + x103;
   const float x187 = -x185*H[8] - x186*H[1];
   const float x188 = -x185*H[9] - x186*H[2];
   const float x189 = -x185*H[10] - x186*H[3];
   const float x190 = -x185*H[12] - x186*H[5];
   const float x191 = -x185*H[13] - x186*H[6];
   const float x192 = -x185*H[14] - x186*H[7];
   const float x193 = -x185*H[11] - x186*H[4] + 1;
   const float x194 = x104 + x105;
   const float x195 = x106 + x107;
   const float x196 = -x194*H[8] - x195*H[1];
   const float x197 = -x194*H[9] - x195*H[2];
   const float x198 = -x194*H[10] - x195*H[3];
   const float x199 = -x194*H[11] - x195*H[4];
   const float x200 = -x194*H[13] - x195*H[6];
   const float x201 = -x194*H[14] - x195*H[7];
   const float x202 = -x194*H[12] - x195*H[5] + 1;
   const float x203 = x108 + x109;
   const float x204 = x110 + x111;
   const float x205 = -x203*H[8] - x204*H[1];
   const float x206 = -x203*H[9] - x204*H[2];
   const float x207 = -x203*H[10] - x204*H[3];
   const float x208 = -x203*H[11] - x204*H[4];
   const float x209 = -x203*H[12] - x204*H[5];
   const float x210 = -x203*H[14] - x204*H[7];
   const float x211 = -x203*H[13] - x204*H[6] + 1;
   const float x212 = x112 + x113;
   const float x213 = x114 + x115;
   const float x214 = -x212*H[8] - x213*H[1];
   const float x215 = -x212*H[9] - x213*H[2];
   const float x216 = -x212*H[10] - x213*H[3];
   const float x217 = -x212*H[11] - x213*H[4];
   const float x218 = -x212*H[12] - x213*H[5];
   const float x219 = -x212*H[13] - x213*H[6];
   const float x220 = -x212*H[14] - x213*H[7] + 1;
   const float x221 = x117 + x119;
   const float x222 = x120 + x121;
   const float x223 = x123 + x125;
   const float x224 = x126 + x127;
   const float x225 = x129 + x131;
   const float x226 = x132 + x133;
   const float x227 = x135 + x137 + x138;
   const float x228 = x139 + x140 + x141;
   xnew[0] = x55*(-x58 - x61 - x63) + x6*(-x48 - x52 - x54) + x64*(-x67 - x70 - x72) + x[0];
   xnew[1] = x55*(-x76 - x77 - x78) + x6*(-x73 - x74 - x75) + x64*(-x79 - x80 - x81) + x[1];
   xnew[2] = x55*(-x85 - x86 - x87) + x6*(-x82 - x83 - x84) + x64*(-x88 - x89 - x90) + x[2];
   xnew[3] = x55*(-x94 - x95 - x96) + x6*(-x91 - x92 - x93) + x64*(-x97 - x98 - x99) + x[3];
   xnew[4] = x55*(-x102 - x103) + x6*(-x100 - x101) + x64*(-x22*x69 - x32*x66) + x[4];
   xnew[5] = x55*(-x106 - x107) + x6*(-x104 - x105) + x64*(-x23*x69 - x33*x66) + x[5];
   xnew[6] = x55*(-x110 - x111) + x6*(-x108 - x109) + x64*(-x24*x69 - x34*x66) + x[6];
   xnew[7] = x55*(-x114 - x115) + x6*(-x112 - x113) + x64*(-x25*x69 - x35*x66) + x[7];
   xnew[8] = x55*(-x120 - x121) + x6*(-x117 - x119) + x64*(-x116*x69 - x118*x66) + x[8];
   xnew[9] = x55*(-x126 - x127) + x6*(-x123 - x125) + x64*(-x122*x69 - x124*x66) + x[9];
   xnew[10] = x55*(-x132 - x133) + x6*(-x129 - x131) + x64*(-x128*x69 - x130*x66) + x[10];
   xnew[11] = x[11];
   xnew[12] = x[12];
   xnew[13] = x[13];
   xnew[14] = x55*(-x139 - x140 - x141) + x6*(-x135 - x137 - x138) + x64*(-x142 - x143 - x144) + x[14];
   Pnew[0] = -x145*P[11] + x148*P[1] + x149*P[2] + x150*P[3] + x151*P[4] + x152*P[5] + x153*P[6] + x154*P[7] + P[0];
   Pnew[1] = -x145*P[22] + x148*P[12] + x149*P[13] + x150*P[14] + x151*P[15] + x152*P[16] + x153*P[17] + x154*P[18] + P[1];
   Pnew[2] = -x145*P[32] + x148*P[13] + x149*P[23] + x150*P[24] + x151*P[25] + x152*P[26] + x153*P[27] + x154*P[28] + P[2];
   Pnew[3] = -x145*P[41] + x148*P[14] + x149*P[24] + x150*P[33] + x151*P[34] + x152*P[35] + x153*P[36] + x154*P[37] + P[3];
   Pnew[4] = x148*P[15] + x149*P[25] + x150*P[34] + x151*P[42] + x152*P[43] + x153*P[44] + x154*P[45] + P[4];
   Pnew[5] = x148*P[16] + x149*P[26] + x150*P[35] + x151*P[43] + x152*P[49] + x153*P[50] + x154*P[51] + P[5];
   Pnew[6] = x148*P[17] + x149*P[27] + x150*P[36] + x151*P[44] + x152*P[50] + x153*P[55] + x154*P[56] + P[6];
   Pnew[7] = x148*P[18] + x149*P[28] + x150*P[37] + x151*P[45] + x152*P[51] + x153*P[56] + x154*P[60] + P[7];
   Pnew[8] = x148*P[19] + x149*P[29] + x150*P[38] + x151*P[46] + x152*P[52] + x153*P[57] + x154*P[61] + P[8];
   Pnew[9] = x148*P[20] + x149*P[30] + x150*P[39] + x151*P[47] + x152*P[53] + x153*P[58] + x154*P[62] + P[9];
   Pnew[10] = x148*P[21] + x149*P[31] + x150*P[40] + x151*P[48] + x152*P[54] + x153*P[59] + x154*P[63] + P[10];
   Pnew[11] = -x145*P[73] + x148*P[22] + x149*P[32] + x150*P[41] + P[11];
   Pnew[12] = -x155*P[22] + x158*P[13] + x159*P[14] + x160*P[15] + x161*P[16] + x162*P[17] + x163*P[18] + x164*P[12];
   Pnew[13] = -x155*P[32] + x158*P[23] + x159*P[24] + x160*P[25] + x161*P[26] + x162*P[27] + x163*P[28] + x164*P[13];
   Pnew[14] = -x155*P[41] + x158*P[24] + x159*P[33] + x160*P[34] + x161*P[35] + x162*P[36] + x163*P[37] + x164*P[14];
   Pnew[15] = x158*P[25] + x159*P[34] + x160*P[42] + x161*P[43] + x162*P[44] + x163*P[45] + x164*P[15];
   Pnew[16] = x158*P[26] + x159*P[35] + x160*P[43] + x161*P[49] + x162*P[50] + x163*P[51] + x164*P[16];
   Pnew[17] = x158*P[27] + x159*P[36] + x160*P[44] + x161*P[50] + x162*P[55] + x163*P[56] + x164*P[17];
   Pnew[18] = x158*P[28] + x159*P[37] + x160*P[45] + x161*P[51] + x162*P[56] + x163*P[60] + x164*P[18];
   Pnew[19] = x158*P[29] + x159*P[38] + x160*P[46] + x161*P[52] + x162*P[57] + x163*P[61] + x164*P[19];
   Pnew[20] = x158*P[30] + x159*P[39] + x160*P[47] + x161*P[53] + x162*P[58] + x163*P[62] + x164*P[20];
   Pnew[21] = x158*P[31] + x159*P[40] + x160*P[48] + x161*P[54] + x162*P[59] + x163*P[63] + x164*P[21];
   Pnew[22] = -x155*P[73] + x158*P[32] + x159*P[41] + x164*P[22];
   Pnew[23] = -x165*P[32] + x168*P[13] + x169*P[24] + x170*P[25] + x171*P[26] + x172*P[27] + x173*P[28] + x174*P[23];
   Pnew[24] = -x165*P[41] + x168*P[14] + x169*P[33] + x170*P[34] + x171*P[35] + x172*P[36] + x173*P[37] + x174*P[24];
   Pnew[25] = x168*P[15] + x169*P[34] + x170*P[42] + x171*P[43] + x172*P[44] + x173*P[45] + x174*P[25];
   Pnew[26] = x168*P[16] + x169*P[35] + x170*P[43] + x171*P[49] + x172*P[50] + x173*P[51] + x174*P[26];
   Pnew[27] = x168*P[17] + x169*P[36] + x170*P[44] + x171*P[50] + x172*P[55] + x173*P[56] + x174*P[27];
   Pnew[28] = x168*P[18] + x169*P[37] + x170*P[45] + x171*P[51] + x172*P[56] + x173*P[60] + x174*P[28];
   Pnew[29] = x168*P[19] + x169*P[38] + x170*P[46] + x171*P[52] + x172*P[57] + x173*P[61] + x174*P[29];
   Pnew[30] = x168*P[20] + x169*P[39] + x170*P[47] + x171*P[53] + x172*P[58] + x173*P[62] + x174*P[30];
   Pnew[31] = x168*P[21] + x169*P[40] + x170*P[48] + x171*P[54] + x172*P[59] + x173*P[63] + x174*P[31];
   Pnew[32] = -x165*P[73] + x168*P[22] + x169*P[41] + x174*P[32];
   Pnew[33] = -x175*P[41] + x178*P[14] + x179*P[24] + x180*P[34] + x181*P[35] + x182*P[36] + x183*P[37] + x184*P[33];
   Pnew[34] = x178*P[15] + x179*P[25] + x180*P[42] + x181*P[43] + x182*P[44] + x183*P[45] + x184*P[34];
   Pnew[35] = x178*P[16] + x179*P[26] + x180*P[43] + x181*P[49] + x182*P[50] + x183*P[51] + x184*P[35];
   Pnew[36] = x178*P[17] + x179*P[27] + x180*P[44] + x181*P[50] + x182*P[55] + x183*P[56] + x184*P[36];
   Pnew[37] = x178*P[18] + x179*P[28] + x180*P[45] + x181*P[51] + x182*P[56] + x183*P[60] + x184*P[37];
   Pnew[38] = x178*P[19] + x179*P[29] + x180*P[46] + x181*P[52] + x182*P[57] + x183*P[61] + x184*P[38];
   Pnew[39] = x178*P[20] + x179*P[30] + x180*P[47] + x181*P[53] + x182*P[58] + x183*P[62] + x184*P[39];
   Pnew[40] = x178*P[21] + x179*P[31] + x180*P[48] + x181*P[54] + x182*P[59] + x183*P[63] + x184*P[40];
   Pnew[41] = -x175*P[73] + x178*P[22] + x179*P[32] + x184*P[41];
   Pnew[42] = x187*P[15] + x188*P[25] + x189*P[34] + x190*P[43] + x191*P[44] + x192*P[45] + x193*P[42];
   Pnew[43] = x187*P[16] + x188*P[26] + x189*P[35] + x190*P[49] + x191*P[50] + x192*P[51] + x193*P[43];
   Pnew[44] = x187*P[17] + x188*P[27] + x189*P[36] + x190*P[50] + x191*P[55] + x192*P[56] + x193*P[44];
   Pnew[45] = x187*P[18] + x188*P[28] + x189*P[37] + x190*P[51] + x191*P[56] + x192*P[60] + x193*P[45];
   Pnew[46] = x187*P[19] + x188*P[29] + x189*P[38] + x190*P[52] + x191*P[57] + x192*P[61] + x193*P[46];
   Pnew[47] = x187*P[20] + x188*P[30] + x189*P[39] + x190*P[53] + x191*P[58] + x192*P[62] + x193*P[47];
   Pnew[48] = x187*P[21] + x188*P[31] + x189*P[40] + x190*P[54] + x191*P[59] + x192*P[63] + x193*P[48];
   Pnew[49] = x196*P[16] + x197*P[26] + x198*P[35] + x199*P[43] + x200*P[50] + x201*P[51] + x202*P[49];
   Pnew[50] = x196*P[17] + x197*P[27] + x198*P[36] + x199*P[44] + x200*P[55] + x201*P[56] + x202*P[50];
   Pnew[51] = x196*P[18] + x197*P[28] + x198*P[37] + x199*P[45] + x200*P[56] + x201*P[60] + x202*P[51];
   Pnew[52] = x196*P[19] + x197*P[29] + x198*P[38] + x199*P[46] + x200*P[57] + x201*P[61] + x202*P[52];
   Pnew[53] = x196*P[20] + x197*P[30] + x198*P[39] + x199*P[47] + x200*P[58] + x201*P[62] + x202*P[53];
   Pnew[54] = x196*P[21] + x197*P[31] + x198*P[40] + x199*P[48] + x200*P[59] + x201*P[63] + x202*P[54];
   Pnew[55] = x205*P[17] + x206*P[27] + x207*P[36] + x208*P[44] + x209*P[50] + x210*P[56] + x211*P[55];
   Pnew[56] = x205*P[18] + x206*P[28] + x207*P[37] + x208*P[45] + x209*P[51] + x210*P[60] + x211*P[56];
   Pnew[57] = x205*P[19] + x206*P[29] + x207*P[38] + x208*P[46] + x209*P[52] + x210*P[61] + x211*P[57];
   Pnew[58] = x205*P[20] + x206*P[30] + x207*P[39] + x208*P[47] + x209*P[53] + x210*P[62] + x211*P[58];
   Pnew[59] = x205*P[21] + x206*P[31] + x207*P[40] + x208*P[48] + x209*P[54] + x210*P[63] + x211*P[59];
   Pnew[60] = x214*P[18] + x215*P[28] + x216*P[37] + x217*P[45] + x218*P[51] + x219*P[56] + x220*P[60];
   Pnew[61] = x214*P[19] + x215*P[29] + x216*P[38] + x217*P[46] + x218*P[52] + x219*P[57] + x220*P[61];
   Pnew[62] = x214*P[20] + x215*P[30] + x216*P[39] + x217*P[47] + x218*P[53] + x219*P[58] + x220*P[62];
   Pnew[63] = x214*P[21] + x215*P[31] + x216*P[40] + x217*P[48] + x218*P[54] + x219*P[59] + x220*P[63];
   Pnew[64] = (-x221*H[8] - x222*H[1])*P[19] + (-x221*H[9] - x222*H[2])*P[29] + (-x221*H[10] - x222*H[3])*P[38] + (-x221*H[11] - x222*H[4])*P[46] + (-x221*H[12] - x222*H[5])*P[52] + (-x221*H[13] - x222*H[6])*P[57] + (-x221*H[14] - x222*H[7])*P[61] + P[64];
   Pnew[65] = P[65];
   Pnew[66] = (-x223*H[8] - x224*H[1])*P[20] + (-x223*H[9] - x224*H[2])*P[30] + (-x223*H[10] - x224*H[3])*P[39] + (-x223*H[11] - x224*H[4])*P[47] + (-x223*H[12] - x224*H[5])*P[53] + (-x223*H[13] - x224*H[6])*P[58] + (-x223*H[14] - x224*H[7])*P[62] + P[66];
   Pnew[67] = P[67];
   Pnew[68] = (-x225*H[8] - x226*H[1])*P[21] + (-x225*H[9] - x226*H[2])*P[31] + (-x225*H[10] - x226*H[3])*P[40] + (-x225*H[11] - x226*H[4])*P[48] + (-x225*H[12] - x226*H[5])*P[54] + (-x225*H[13] - x226*H[6])*P[59] + (-x225*H[14] - x226*H[7])*P[63] + P[68];
   Pnew[69] = P[69];
   Pnew[70] = P[70];
   Pnew[71] = P[71];
   Pnew[72] = P[72];
   Pnew[73] = (-x227*H[8] - x228*H[1])*P[22] + (-x227*H[9] - x228*H[2])*P[32] + (-x227*H[10] - x228*H[3])*P[41] + (-(x142 + x143 + x144)*H[15] + 1)*P[73];
}

void covariance_init(float *Pnew) {
	float P0[NUMX] = {1.0e7f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		              1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
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
	Pnew[49] = P0[5];
	Pnew[50] = 0;
	Pnew[51] = 0;
	Pnew[52] = 0;
	Pnew[53] = 0;
	Pnew[54] = 0;
	Pnew[55] = P0[6];
	Pnew[56] = 0;
	Pnew[57] = 0;
	Pnew[58] = 0;
	Pnew[59] = 0;
	Pnew[60] = P0[7];
	Pnew[61] = 0;
	Pnew[62] = 0;
	Pnew[63] = 0;
	Pnew[64] = P0[8];
	Pnew[65] = 0;
	Pnew[66] = P0[9];
	Pnew[67] = 0;
	Pnew[68] = P0[10];
	Pnew[69] = 0;
	Pnew[70] = P0[11];
	Pnew[71] = P0[12];
	Pnew[72] = P0[13];
	Pnew[73] = P0[14];
}

/**
 * @}
 * @}
 */
