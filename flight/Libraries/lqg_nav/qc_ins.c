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
#include "stdint.h"

#include "stdio.h"

// Local definitions
#define NUMX 15  // number of state variables
#define NUMF 56  // number of non-zero variables stored in the F matrix
#define NUMH 27  // number of non-zero variables stored in the H matrix
#define NUMP 74  // number of non-zero variables in the upper triangular of covariance
#define NUMO 7   // number of total outputs

// Private methods
void covariance_init(double *Pnew);
void state_prediction(double *x, const float *u, double Ts, double *param, double *xnew);
void covariance_prediction(double *P, double Ts, double *F, double *Q, double *Pnew);
void linearize_FH(double *x, double *param, double *F, double *H);
//void baro_correction(double *x, double *P, double baro, double *H, double *R, double *xnew, double *Pnew);
void gyro_correction(double *x, double *P, const float *gyro, double *H, double *R, double *xnew, double *Pnew);
//void accel_correction(double *x, double *P, double *accel, double *H, double *R, double *xnew, double *Pnew);
void update_state(double *x, double *P, double *xnew, double *Pnew);
void normalize_state(double *x);

enum qcins_state_magic {
  QCINS_STATE_MAGIC = 0xe7b5500e, // echo "qc_ins.c" | md5
};

struct qcins_state {

	double Ts;      // The kalman filter time step

	struct {
		double g;
		double beta1;
		double tau;
		double mu;
	} __attribute__((__packed__)) params;

	double x[NUMX]; // buffer to store the current state
	double P[NUMP]; // buffer to store the covariance
	double F[NUMF]; // buffer to store the state dynamics derivatives
	double H[NUMH]; // buffer to store the output derivatives
	double Q[NUMX]; // process covariance
	double R[NUMO]; // measurement noises

	double xnew[NUMX];
	double Pnew[NUMP];

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
		qcins_state->Q[i] = 1e-5;
	qcins_state->Q[NUMX-1] = 1;
	qcins_state->Q[NUMX-2] = 1;
	qcins_state->Q[NUMX-3] = 1;
	qcins_state->Q[NUMX-4] = 1;

	// Store the observation noises
	for (uint32_t i = 0; i < NUMX; i++)
		qcins_state->R[i] = 1000.0f;

	return true;
}

bool qcins_alloc(uintptr_t *qcins_handle)
{
	struct qcins_state *qcins_state = (struct qcins_state *) PIOS_malloc(sizeof(struct qcins_state));
	if (qcins_state == NULL)
		return false;

	qcins_state->magic = QCINS_STATE_MAGIC;

	qcins_state->params.g = 9.81f;
	qcins_state->params.beta1 = 1000.0f;
	qcins_state->params.tau = 0.015f;
	qcins_state->params.mu = 1;

	qcins_init((uintptr_t) qcins_state);

	(*qcins_handle) = (uintptr_t) qcins_state;

	return true;
}

bool qcins_predict(uintptr_t qcins_handle, const float roll, const float pitch, const float yaw, const float throttle, double Ts)
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

bool qcins_correct_accel_gyro(uintptr_t qcins_handle, const float gyros[3], const float accels[3])
{
	struct qcins_state *qcins_state = (struct qcins_state *) qcins_handle;
	if (!qcins_validate(qcins_state))
		return false;

	gyro_correction(qcins_state->x, qcins_state->P, gyros, qcins_state->H, qcins_state->R, qcins_state->xnew, qcins_state->Pnew);
	/*if (0)
		accel_correction(qcins_state->xnew, qcins_state->Pnew, accels, qcins_state->H, qcins_state->R, qcins_state->x, qcins_state->P);
	else*/
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

void normalize_state(double *x)
{
	double qmag = sqrtf(x[4]*x[4] + x[5]*x[5] + x[6]*x[6] + x[7]*x[7]);
	x[4] /= qmag;
	x[5] /= qmag;
	x[6] /= qmag;
	x[7] /= qmag;
}

void update_state(double *x, double *P, double *xnew, double *Pnew)
{
	for (int i = 0; i < NUMX; i++)
		xnew[i] = x[i];
	for (int i = 0; i < NUMP; i++)
		Pnew[i] = P[i];
}

void state_prediction(double *x, const float *u, double Ts, double *param, double *xnew) {
   const double x0 = 2*param[0]*x[14];
   const double x1 = x[4]*x[6];
   const double x2 = x[5]*x[7];
   const double x3 = 2*x[4]*x[7];
   const double x4 = x3 - 2*x[5]*x[6];
   const double x5 = 2*x[3];
   const double x6 = x[4]*x[5];
   const double x7 = x[6]*x[7];
   const double x8 = x6 + x7;
   const double x9 = (x[4]*x[4]);
   const double x10 = (x[5]*x[5]);
   const double x11 = -x10 + x9;
   const double x12 = (x[6]*x[6]);
   const double x13 = (x[7]*x[7]);
   const double x14 = -x13;
   const double x15 = x11 + x12 + x14;
   const double x16 = x15*x[2] - x4*x[1] + x5*x8;
   const double x17 = x16*param[3];
   const double x18 = -x12;
   const double x19 = x10 + x14 + x18 + x9;
   const double x20 = x3 + 2*x[5]*x[6];
   const double x21 = x1 - x2;
   const double x22 = x19*x[1] + x20*x[2] - x21*x5;
   const double x23 = x22*param[3];
   const double x24 = 2*param[3];
   const double x25 = (1.0F/2.0F)*x[5];
   const double x26 = (1.0F/2.0F)*x[6];
   const double x27 = (1.0F/2.0F)*x[7];
   const double x28 = (1.0F/2.0F)*x[4];
   const double x29 = Ts*param[1];
   const double x30 = Ts/param[2];
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

void linearize_FH(double *x, double *param, double *F, double *H) {
   const double x0 = 2*x[4];
   const double x1 = x0*x[7];
   const double x2 = 2*x[5];
   const double x3 = x2*x[6];
   const double x4 = -x1 + x3;
   const double x5 = (x[4]*x[4]);
   const double x6 = (x[5]*x[5]);
   const double x7 = (x[6]*x[6]);
   const double x8 = -x7;
   const double x9 = (x[7]*x[7]);
   const double x10 = -x9;
   const double x11 = x10 + x5 + x6 + x8;
   const double x12 = x4*param[3];
   const double x13 = x5 - x6;
   const double x14 = x10 + x13 + x7;
   const double x15 = x1 + x3;
   const double x16 = x15*param[3];
   const double x17 = -x11*x16 - x12*x14;
   const double x18 = x0*x[5];
   const double x19 = 2*x[6];
   const double x20 = x19*x[7];
   const double x21 = x18 + x20;
   const double x22 = x21*param[3];
   const double x23 = x0*x[6];
   const double x24 = x2*x[7];
   const double x25 = -x23 + x24;
   const double x26 = x25*param[3];
   const double x27 = -x11*x26 - x22*x4;
   const double x28 = 2*x[7];
   const double x29 = x0*x[2] + x2*x[3] - x28*x[1];
   const double x30 = x0*x[1];
   const double x31 = x28*x[2];
   const double x32 = x19*x[3];
   const double x33 = x30 + x31 - x32;
   const double x34 = x11*param[3];
   const double x35 = param[0]*x[14];
   const double x36 = x19*x35;
   const double x37 = (x14*x[2] + x21*x[3] + x4*x[1])*param[3];
   const double x38 = x28*x37;
   const double x39 = (x11*x[1] + x15*x[2] + x25*x[3])*param[3];
   const double x40 = x0*x39;
   const double x41 = -x36 + x38 - x40;
   const double x42 = x19*x[1];
   const double x43 = x2*x[2];
   const double x44 = x0*x[3];
   const double x45 = x42 - x43 + x44;
   const double x46 = x19*x[2] + x2*x[1] + x28*x[3];
   const double x47 = -x19*x37 - x2*x39 - x28*x35;
   const double x48 = -x42 + x43 - x44;
   const double x49 = x0*x35;
   const double x50 = x2*x37;
   const double x51 = x19*x39;
   const double x52 = -x49 - x50 + x51;
   const double x53 = x2*x35;
   const double x54 = -x30 - x31 + x32;
   const double x55 = x0*x37;
   const double x56 = x28*x39;
   const double x57 = -x14*x22 - x15*x26;
   const double x58 = x14*param[3];
   const double x59 = x53 - x55 - x56;
   const double x60 = (1.0F/2.0F)*x[8];
   const double x61 = -x60;
   const double x62 = (1.0F/2.0F)*x[9];
   const double x63 = -x62;
   const double x64 = (1.0F/2.0F)*x[10];
   const double x65 = -x64;
   const double x66 = (1.0F/2.0F)*x[5];
   const double x67 = -x66;
   const double x68 = (1.0F/2.0F)*x[6];
   const double x69 = -x68;
   const double x70 = (1.0F/2.0F)*x[7];
   const double x71 = -x70;
   const double x72 = (1.0F/2.0F)*x[4];
   const double x73 = -1/param[2];
   const double x74 = -x46*param[3];
   const double x75 = -x29*param[3];
   const double x76 = -x28;
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
void covariance_prediction(double *P, double Ts, double *F, double *Q, double *Pnew) {
   const double x0 = Ts*F[0];
   const double x1 = x0*P[33] + P[3];
   const double x2 = Ts*F[2];
   const double x3 = x0*P[24] + P[2];
   const double x4 = Ts*F[3];
   const double x5 = Ts*F[4];
   const double x6 = x0*P[34] + P[4];
   const double x7 = Ts*F[5];
   const double x8 = x0*P[35] + P[5];
   const double x9 = Ts*F[6];
   const double x10 = x0*P[36] + P[6];
   const double x11 = Ts*F[7];
   const double x12 = x0*P[37] + P[7];
   const double x13 = Ts*F[8];
   const double x14 = x0*P[41] + P[11];
   const double x15 = Ts*F[1] + 1;
   const double x16 = x0*P[14] + P[1];
   const double x17 = Ts*F[9];
   const double x18 = Ts*F[11];
   const double x19 = Ts*F[12];
   const double x20 = Ts*F[13];
   const double x21 = Ts*F[14];
   const double x22 = Ts*F[15];
   const double x23 = Ts*F[16];
   const double x24 = Ts*F[10] + 1;
   const double x25 = Ts*F[17];
   const double x26 = Ts*F[18];
   const double x27 = Ts*F[20];
   const double x28 = Ts*F[21];
   const double x29 = Ts*F[22];
   const double x30 = Ts*F[23];
   const double x31 = Ts*F[24];
   const double x32 = Ts*F[19] + 1;
   const double x33 = Ts*F[25];
   const double x34 = Ts*F[26];
   const double x35 = Ts*F[27];
   const double x36 = Ts*F[28];
   const double x37 = x0*P[38] + P[8];
   const double x38 = Ts*F[29];
   const double x39 = x0*P[39] + P[9];
   const double x40 = Ts*F[30];
   const double x41 = x0*P[40] + P[10];
   const double x42 = Ts*F[31];
   const double x43 = Ts*F[32];
   const double x44 = Ts*F[33];
   const double x45 = Ts*F[34];
   const double x46 = Ts*F[35];
   const double x47 = Ts*F[36];
   const double x48 = Ts*F[37];
   const double x49 = Ts*F[38];
   const double x50 = Ts*F[39];
   const double x51 = Ts*F[40];
   const double x52 = Ts*F[41];
   const double x53 = Ts*F[42];
   const double x54 = Ts*F[43];
   const double x55 = Ts*F[44];
   const double x56 = Ts*F[45];
   const double x57 = Ts*F[46];
   const double x58 = Ts*F[47];
   const double x59 = Ts*F[48];
   const double x60 = Ts*F[55] + 1;
   const double x61 = x13*P[73] + x15*P[22] + x2*P[32] + x4*P[41];
   const double x62 = x11*P[45] + x15*P[15] + x2*P[25] + x4*P[34] + x5*P[42] + x7*P[43] + x9*P[44];
   const double x63 = x11*P[51] + x15*P[16] + x2*P[26] + x4*P[35] + x5*P[43] + x7*P[49] + x9*P[50];
   const double x64 = x11*P[56] + x15*P[17] + x2*P[27] + x4*P[36] + x5*P[44] + x7*P[50] + x9*P[55];
   const double x65 = x11*P[60] + x15*P[18] + x2*P[28] + x4*P[37] + x5*P[45] + x7*P[51] + x9*P[56];
   const double x66 = x11*P[28] + x13*P[32] + x15*P[13] + x2*P[23] + x4*P[24] + x5*P[25] + x7*P[26] + x9*P[27];
   const double x67 = x11*P[37] + x13*P[41] + x15*P[14] + x2*P[24] + x4*P[33] + x5*P[34] + x7*P[35] + x9*P[36];
   const double x68 = x11*P[18] + x13*P[22] + x15*P[12] + x2*P[13] + x4*P[14] + x5*P[15] + x7*P[16] + x9*P[17];
   const double x69 = x11*P[61] + x15*P[19] + x2*P[29] + x4*P[38] + x5*P[46] + x7*P[52] + x9*P[57];
   const double x70 = x11*P[62] + x15*P[20] + x2*P[30] + x4*P[39] + x5*P[47] + x7*P[53] + x9*P[58];
   const double x71 = x11*P[63] + x15*P[21] + x2*P[31] + x4*P[40] + x5*P[48] + x7*P[54] + x9*P[59];
   const double x72 = x17*P[22] + x18*P[41] + x23*P[73] + x24*P[32];
   const double x73 = x17*P[15] + x18*P[34] + x19*P[42] + x20*P[43] + x21*P[44] + x22*P[45] + x24*P[25];
   const double x74 = x17*P[16] + x18*P[35] + x19*P[43] + x20*P[49] + x21*P[50] + x22*P[51] + x24*P[26];
   const double x75 = x17*P[17] + x18*P[36] + x19*P[44] + x20*P[50] + x21*P[55] + x22*P[56] + x24*P[27];
   const double x76 = x17*P[18] + x18*P[37] + x19*P[45] + x20*P[51] + x21*P[56] + x22*P[60] + x24*P[28];
   const double x77 = x17*P[12] + x18*P[14] + x19*P[15] + x20*P[16] + x21*P[17] + x22*P[18] + x23*P[22] + x24*P[13];
   const double x78 = x17*P[14] + x18*P[33] + x19*P[34] + x20*P[35] + x21*P[36] + x22*P[37] + x23*P[41] + x24*P[24];
   const double x79 = x17*P[13] + x18*P[24] + x19*P[25] + x20*P[26] + x21*P[27] + x22*P[28] + x23*P[32] + x24*P[23];
   const double x80 = x17*P[19] + x18*P[38] + x19*P[46] + x20*P[52] + x21*P[57] + x22*P[61] + x24*P[29];
   const double x81 = x17*P[20] + x18*P[39] + x19*P[47] + x20*P[53] + x21*P[58] + x22*P[62] + x24*P[30];
   const double x82 = x17*P[21] + x18*P[40] + x19*P[48] + x20*P[54] + x21*P[59] + x22*P[63] + x24*P[31];
   const double x83 = x25*P[22] + x26*P[32] + x31*P[73] + x32*P[41];
   const double x84 = x25*P[15] + x26*P[25] + x27*P[42] + x28*P[43] + x29*P[44] + x30*P[45] + x32*P[34];
   const double x85 = x25*P[16] + x26*P[26] + x27*P[43] + x28*P[49] + x29*P[50] + x30*P[51] + x32*P[35];
   const double x86 = x25*P[17] + x26*P[27] + x27*P[44] + x28*P[50] + x29*P[55] + x30*P[56] + x32*P[36];
   const double x87 = x25*P[18] + x26*P[28] + x27*P[45] + x28*P[51] + x29*P[56] + x30*P[60] + x32*P[37];
   const double x88 = x25*P[19] + x26*P[29] + x27*P[46] + x28*P[52] + x29*P[57] + x30*P[61] + x32*P[38];
   const double x89 = x25*P[20] + x26*P[30] + x27*P[47] + x28*P[53] + x29*P[58] + x30*P[62] + x32*P[39];
   const double x90 = x25*P[21] + x26*P[31] + x27*P[48] + x28*P[54] + x29*P[59] + x30*P[63] + x32*P[40];
   const double x91 = x33*P[43] + x34*P[44] + x35*P[45] + x36*P[46] + x38*P[47] + x40*P[48] + P[42];
   const double x92 = x33*P[52] + x34*P[57] + x35*P[61] + x36*P[64] + P[46];
   const double x93 = x33*P[53] + x34*P[58] + x35*P[62] + x38*P[66] + P[47];
   const double x94 = x33*P[54] + x34*P[59] + x35*P[63] + x40*P[68] + P[48];
   const double x95 = x33*P[49] + x34*P[50] + x35*P[51] + x36*P[52] + x38*P[53] + x40*P[54] + P[43];
   const double x96 = x33*P[50] + x34*P[55] + x35*P[56] + x36*P[57] + x38*P[58] + x40*P[59] + P[44];
   const double x97 = x33*P[51] + x34*P[56] + x35*P[60] + x36*P[61] + x38*P[62] + x40*P[63] + P[45];
   const double x98 = (Ts*Ts);
   const double x99 = x98*F[49]*P[65];
   const double x100 = x98*F[50]*P[67];
   const double x101 = x98*F[51]*P[69];
   const double x102 = x42*P[43] + x43*P[50] + x44*P[51] + x45*P[52] + x46*P[53] + x47*P[54] + P[49];
   const double x103 = x42*P[46] + x43*P[57] + x44*P[61] + x45*P[64] + P[52];
   const double x104 = x42*P[47] + x43*P[58] + x44*P[62] + x46*P[66] + P[53];
   const double x105 = x42*P[48] + x43*P[59] + x44*P[63] + x47*P[68] + P[54];
   const double x106 = x42*P[42] + x43*P[44] + x44*P[45] + x45*P[46] + x46*P[47] + x47*P[48] + P[43];
   const double x107 = x42*P[44] + x43*P[55] + x44*P[56] + x45*P[57] + x46*P[58] + x47*P[59] + P[50];
   const double x108 = x42*P[45] + x43*P[56] + x44*P[60] + x45*P[61] + x46*P[62] + x47*P[63] + P[51];
   const double x109 = x48*P[44] + x49*P[50] + x50*P[56] + x51*P[57] + x52*P[58] + x53*P[59] + P[55];
   const double x110 = x48*P[46] + x49*P[52] + x50*P[61] + x51*P[64] + P[57];
   const double x111 = x48*P[47] + x49*P[53] + x50*P[62] + x52*P[66] + P[58];
   const double x112 = x48*P[48] + x49*P[54] + x50*P[63] + x53*P[68] + P[59];
   const double x113 = x48*P[42] + x49*P[43] + x50*P[45] + x51*P[46] + x52*P[47] + x53*P[48] + P[44];
   const double x114 = x48*P[43] + x49*P[49] + x50*P[51] + x51*P[52] + x52*P[53] + x53*P[54] + P[50];
   const double x115 = x48*P[45] + x49*P[51] + x50*P[60] + x51*P[61] + x52*P[62] + x53*P[63] + P[56];
   const double x116 = x54*P[46] + x55*P[52] + x56*P[57] + x57*P[64] + P[61];
   const double x117 = x54*P[47] + x55*P[53] + x56*P[58] + x58*P[66] + P[62];
   const double x118 = x54*P[48] + x55*P[54] + x56*P[59] + x59*P[68] + P[63];
   const double x119 = Ts*F[49];
   const double x120 = x119*P[70] + P[65];
   const double x121 = Ts*F[52] + 1;
   const double x122 = Ts*F[50];
   const double x123 = x122*P[71] + P[67];
   const double x124 = Ts*F[53] + 1;
   const double x125 = Ts*F[51];
   const double x126 = x125*P[72] + P[69];
   const double x127 = Ts*F[54] + 1;
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

void baro_correction(double *x, double *P, double baro, double *H, double *R, double *xnew, double *Pnew) {
   const double x0 = (H[0]*H[0]);
   const double x1 = x0*P[0];
   const double x2 = 1.0F/(x1 + R[0]);
   const double x3 = x2*(-baro + x[0])*H[0];
   const double x4 = -x1*x2 + 1;
   const double x5 = x0*x2;
   const double x6 = x0*x2*P[1];
   const double x7 = x0*x2*P[2];
   const double x8 = x0*x2*P[3];
   const double x9 = x0*x2*P[4];
   const double x10 = x0*x2*P[5];
   const double x11 = x0*x2*P[6];
   const double x12 = x0*x2*P[7];
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

void gyro_correction(double *x, double *P, const float *gyro, double *H, double *R, double *xnew, double *Pnew) {
   const double x0 = (H[16]*H[16]);
   const double x1 = x0*P[64];
   const double x2 = 1.0F/(x1 + R[4]);
   const double x3 = x2*(-gyro[0] + x[8])*H[16];
   const double x4 = (H[17]*H[17]);
   const double x5 = x4*P[66];
   const double x6 = 1.0F/(x5 + R[5]);
   const double x7 = x6*(-gyro[1] + x[9])*H[17];
   const double x8 = (H[18]*H[18]);
   const double x9 = x8*P[68];
   const double x10 = 1.0F/(x9 + R[6]);
   const double x11 = x10*(-gyro[2] + x[10])*H[18];
   const double x12 = x0*x2;
   const double x13 = x4*x6;
   const double x14 = x10*x8;
   const double x15 = x0*x2*P[8];
   const double x16 = x4*x6*P[9];
   const double x17 = x10*x8*P[10];
   const double x18 = x1*x2;
   const double x19 = x5*x6;
   const double x20 = x10*x9;
   const double x21 = x0*x2*P[19];
   const double x22 = x4*x6*P[20];
   const double x23 = x10*x8*P[21];
   const double x24 = x0*x2*P[29];
   const double x25 = x4*x6*P[30];
   const double x26 = x10*x8*P[31];
   const double x27 = x0*x2*P[38];
   const double x28 = x4*x6*P[39];
   const double x29 = x10*x8*P[40];
   const double x30 = x0*x2*P[46];
   const double x31 = x4*x6*P[47];
   const double x32 = x10*x8*P[48];
   const double x33 = x0*x2*P[52];
   const double x34 = x4*x6*P[53];
   const double x35 = x10*x8*P[54];
   const double x36 = -x18 + 1;
   const double x37 = -x19 + 1;
   const double x38 = -x20 + 1;
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

void covariance_init(double *Pnew) {
	float P0[NUMX] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
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