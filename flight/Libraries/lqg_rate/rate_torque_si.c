/**
 ******************************************************************************
 * @addtogroup TauLabsModules Tau Labs Modules
 * @{
 * @addtogroup Rate Torque Linear Quadratic Gaussian controller
 * @{
 *
 * @file       rate_torque_si.c
 * @author     Tau Labs, http://taulabs.org, Copyright (C) 2016
 * @brief      System identication of parameters of this model
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

#include "math.h"
#include "stdint.h"
#include "rate_torque_si.h"

 /**
 * Prediction step for EKF on control inputs to quad that
 * learns the system properties
 * @param X the current state estimate which is updated in place
 * @param P the current covariance matrix, updated in place
 * @param[in] the current control inputs (roll, pitch, yaw)
 * @param[in] the gyro measurements
 */
/**
 * Prediction step for EKF on control inputs to quad that
 * learns the system properties
 * @param X the current state estimate which is updated in place
 * @param P the current covariance matrix, updated in place
 * @param[in] the current control inputs (roll, pitch, yaw)
 * @param[in] the gyro measurements
 */
void rtsi_predict(float X[RTSI_NUMX], float P[RTSI_NUMP], const float u_in[3], const float gyro[3], const float dT_s)
{

	const float Ts = dT_s;
	const float Tsq = Ts * Ts;
	//const float Tsq3 = Tsq * Ts;
	//const float Tsq4 = Tsq * Tsq;

	// for convenience and clarity code below uses the named versions of
	// the state variables
	float w1 = X[0];           // roll rate estimate
	float w2 = X[1];           // pitch rate estimate
	float w3 = X[2];           // yaw rate estimate
	float u1 = X[3];           // scaled roll torque 
	float u2 = X[4];           // scaled pitch torque
	float u3 = X[5];           // scaled yaw torque
	const float e_b1 = expf(X[6]);   // roll torque scale
	const float b1 = X[6];
	const float e_b2 = expf(X[7]);   // pitch torque scale
	const float b2 = X[7];
	const float e_b3 = expf(X[8]);   // yaw torque scale
	const float b3 = X[8];
	const float e_tau = expf(X[9]); // time response of the motors
	const float tau = X[9];
	const float e_tau_y = expf(X[10]); // time response of the motors
	const float tau_y = X[10];
	const float bias1 = X[11];        // bias in the roll torque
	const float bias2 = X[12];       // bias in the pitch torque
	const float bias3 = X[13];       // bias in the yaw torque

	// inputs to the system (roll, pitch, yaw)
	const float u1_in = u_in[0];
	const float u2_in = u_in[1];
	const float u3_in = u_in[2];

	// measurements from gyro
	const float gyro_x = gyro[0];
	const float gyro_y = gyro[1];
	const float gyro_z = gyro[2];

	// update named variables because we want to use predicted
	// values below
	w1 = X[0] = w1 - Ts*bias1*e_b1 + Ts*u1*e_b1;
	w2 = X[1] = w2 - Ts*bias2*e_b2 + Ts*u2*e_b2;
	w3 = X[2] = w3 - Ts*bias3*e_b3 + Ts*u3*e_b3;
	u1 = X[3] = (Ts*u1_in)/(Ts + e_tau) + (u1*e_tau)/(Ts + e_tau);
	u2 = X[4] = (Ts*u2_in)/(Ts + e_tau) + (u2*e_tau)/(Ts + e_tau);
	u3 = X[5] = (Ts*u3_in)/(Ts + e_tau_y) + (u3*e_tau_y)/(Ts + e_tau_y);
    // X[6] to X[13] unchanged

	/**** filter parameters ****/
	// core state variables, these were determined from offline analysis and replay of flights
	const float q_w = 1e0f;
	const float q_ud = 1e-5f;
	const float q_bias = 1e-19f;
	const float s_a = 3000.0f;  // expected gyro noise
	// system identification parameters
	const float q_B = 1e-5f;
	const float q_tau = 1e-5f;

	const float Q[RTSI_NUMX] = {q_w, q_w, q_w, q_ud, q_ud, q_ud, q_B, q_B, q_B, q_tau, q_tau, q_bias, q_bias, q_bias};

	float D[RTSI_NUMP];
	for (uint32_t i = 0; i < RTSI_NUMP; i++)
        D[i] = P[i];

    const float e_tau2    = e_tau * e_tau;
    //const float e_tau3    = e_tau * e_tau2;
    //const float e_tau4    = e_tau2 * e_tau2;
    const float Ts_e_tau2 = (Ts + e_tau) * (Ts + e_tau);
    const float Ts_e_tau3 = Ts_e_tau2 * (Ts + e_tau);
    const float Ts_e_tau4 = Ts_e_tau2 * Ts_e_tau2;
    const float e_tau_y2    = e_tau_y * e_tau_y;
    //const float e_tau_y3    = e_tau_y * e_tau_y2;
    //const float e_tau_y4    = e_tau_y2 * e_tau_y2;
    const float Ts_e_tau_y2 = (Ts + e_tau_y) * (Ts + e_tau_y);
    const float Ts_e_tau_y3 = Ts_e_tau_y2 * (Ts + e_tau_y);
    const float Ts_e_tau_y4 = Ts_e_tau_y2 * Ts_e_tau_y2;

	// covariance propagation - D is stored copy of covariance	
	P[0] = D[0] + Q[0] + D[4]*Tsq*(e_b1*e_b1) - 2*D[30]*Tsq*(e_b1*e_b1) + D[33]*Tsq*(e_b1*e_b1) + 2*D[3]*Ts*e_b1 - 2*D[29]*Ts*e_b1 - 2*D[9]*Ts*bias1*e_b1 + 2*D[9]*Ts*u1*e_b1 + D[11]*Tsq*(bias1*bias1)*(e_b1*e_b1) + D[11]*Tsq*(u1*u1)*(e_b1*e_b1) - 2*D[10]*Tsq*bias1*(e_b1*e_b1) + 2*D[31]*Tsq*bias1*(e_b1*e_b1) + 2*D[10]*Tsq*u1*(e_b1*e_b1) - 2*D[31]*Tsq*u1*(e_b1*e_b1) - 2*D[11]*Tsq*bias1*u1*(e_b1*e_b1);
	P[1] = D[1] + Q[1] + D[6]*Tsq*(e_b2*e_b2) - 2*D[35]*Tsq*(e_b2*e_b2) + D[38]*Tsq*(e_b2*e_b2) + 2*D[5]*Ts*e_b2 - 2*D[34]*Ts*e_b2 - 2*D[12]*Ts*bias2*e_b2 + 2*D[12]*Ts*u2*e_b2 + D[14]*Tsq*(bias2*bias2)*(e_b2*e_b2) + D[14]*Tsq*(u2*u2)*(e_b2*e_b2) - 2*D[13]*Tsq*bias2*(e_b2*e_b2) + 2*D[36]*Tsq*bias2*(e_b2*e_b2) + 2*D[13]*Tsq*u2*(e_b2*e_b2) - 2*D[36]*Tsq*u2*(e_b2*e_b2) - 2*D[14]*Tsq*bias2*u2*(e_b2*e_b2);
	P[2] = D[2] + Q[2] + D[8]*Tsq*(e_b3*e_b3) - 2*D[40]*Tsq*(e_b3*e_b3) + D[43]*Tsq*(e_b3*e_b3) + 2*D[7]*Ts*e_b3 - 2*D[39]*Ts*e_b3 - 2*D[15]*Ts*bias3*e_b3 + 2*D[15]*Ts*u3*e_b3 + D[17]*Tsq*(bias3*bias3)*(e_b3*e_b3) + D[17]*Tsq*(u3*u3)*(e_b3*e_b3) - 2*D[16]*Tsq*bias3*(e_b3*e_b3) + 2*D[41]*Tsq*bias3*(e_b3*e_b3) + 2*D[16]*Tsq*u3*(e_b3*e_b3) - 2*D[41]*Tsq*u3*(e_b3*e_b3) - 2*D[17]*Tsq*bias3*u3*(e_b3*e_b3);
	P[3] = (Ts*e_tau*(u1 - u1_in)*(D[18] + D[20]*Ts*e_b1 - D[32]*Ts*e_b1 - D[22]*Ts*e_b1*(bias1 - u1)))/Ts_e_tau2 - (Ts/(Ts + e_tau) - 1)*(D[3] + D[4]*Ts*e_b1 - D[30]*Ts*e_b1 - D[10]*Ts*e_b1*(bias1 - u1));
	P[4] = Q[3] + (e_tau2*(D[4]*Ts + D[4]*e_tau + D[20]*Ts*u1 - D[20]*Ts*u1_in))/Ts_e_tau3 + (Ts*e_tau2*(u1 - u1_in)*(D[20]*Ts + D[20]*e_tau + D[24]*Ts*u1 - D[24]*Ts*u1_in))/Ts_e_tau4;
	P[5] = (Ts*e_tau*(u2 - u2_in)*(D[19] + D[21]*Ts*e_b2 - D[37]*Ts*e_b2 - D[23]*Ts*e_b2*(bias2 - u2)))/Ts_e_tau2 - (Ts/(Ts + e_tau) - 1)*(D[5] + D[6]*Ts*e_b2 - D[35]*Ts*e_b2 - D[13]*Ts*e_b2*(bias2 - u2));
	P[6] = Q[4] + (e_tau2*(D[6]*Ts + D[6]*e_tau + D[21]*Ts*u2 - D[21]*Ts*u2_in))/Ts_e_tau3 + (Ts*e_tau2*(u2 - u2_in)*(D[21]*Ts + D[21]*e_tau + D[24]*Ts*u2 - D[24]*Ts*u2_in))/Ts_e_tau4;
	P[7] = (Ts*e_tau_y*(u3 - u3_in)*(D[25] + D[26]*Ts*e_b3 - D[42]*Ts*e_b3 - D[27]*Ts*e_b3*(bias3 - u3)))/Ts_e_tau_y2 - (Ts/(Ts + e_tau_y) - 1)*(D[7] + D[8]*Ts*e_b3 - D[40]*Ts*e_b3 - D[16]*Ts*e_b3*(bias3 - u3));
	P[8] = Q[5] + (e_tau_y2*(D[8]*Ts + D[8]*e_tau_y + D[26]*Ts*u3 - D[26]*Ts*u3_in))/Ts_e_tau_y3 + (Ts*e_tau_y2*(u3 - u3_in)*(D[26]*Ts + D[26]*e_tau_y + D[28]*Ts*u3 - D[28]*Ts*u3_in))/Ts_e_tau_y4;
	P[9] = D[9] + D[10]*Ts*e_b1 - D[31]*Ts*e_b1 - D[11]*Ts*e_b1*(bias1 - u1);
	P[10] = (e_tau*(D[10]*Ts + D[10]*e_tau + D[22]*Ts*u1 - D[22]*Ts*u1_in))/Ts_e_tau2;
	P[11] = D[11] + Q[6];
	P[12] = D[12] + D[13]*Ts*e_b2 - D[36]*Ts*e_b2 - D[14]*Ts*e_b2*(bias2 - u2);
	P[13] = (e_tau*(D[13]*Ts + D[13]*e_tau + D[23]*Ts*u2 - D[23]*Ts*u2_in))/Ts_e_tau2;
	P[14] = D[14] + Q[7];
	P[15] = D[15] + D[16]*Ts*e_b3 - D[41]*Ts*e_b3 - D[17]*Ts*e_b3*(bias3 - u3);
	P[16] = (e_tau_y*(D[16]*Ts + D[16]*e_tau_y + D[27]*Ts*u3 - D[27]*Ts*u3_in))/Ts_e_tau_y2;
	P[17] = D[17] + Q[8];
	P[18] = D[18] + D[20]*Ts*e_b1 - D[32]*Ts*e_b1 - D[22]*Ts*e_b1*(bias1 - u1);
	P[19] = D[19] + D[21]*Ts*e_b2 - D[37]*Ts*e_b2 - D[23]*Ts*e_b2*(bias2 - u2);
	P[20] = (e_tau*(D[20]*Ts + D[20]*e_tau + D[24]*Ts*u1 - D[24]*Ts*u1_in))/Ts_e_tau2;
	P[21] = (e_tau*(D[21]*Ts + D[21]*e_tau + D[24]*Ts*u2 - D[24]*Ts*u2_in))/Ts_e_tau2;
	P[22] = D[22];
	P[23] = D[23];
	P[24] = D[24] + Q[9];
	P[25] = D[25] + D[26]*Ts*e_b3 - D[42]*Ts*e_b3 - D[27]*Ts*e_b3*(bias3 - u3);
	P[26] = (e_tau_y*(D[26]*Ts + D[26]*e_tau_y + D[28]*Ts*u3 - D[28]*Ts*u3_in))/Ts_e_tau_y2;
	P[27] = D[27];
	P[28] = D[28] + Q[10];
	P[29] = D[29] + D[30]*Ts*e_b1 - D[33]*Ts*e_b1 - D[31]*Ts*e_b1*(bias1 - u1);
	P[30] = (e_tau*(D[30]*Ts + D[30]*e_tau + D[32]*Ts*u1 - D[32]*Ts*u1_in))/Ts_e_tau2;
	P[31] = D[31];
	P[32] = D[32];
	P[33] = D[33] + Q[11];
	P[34] = D[34] + D[35]*Ts*e_b2 - D[38]*Ts*e_b2 - D[36]*Ts*e_b2*(bias2 - u2);
	P[35] = (e_tau*(D[35]*Ts + D[35]*e_tau + D[37]*Ts*u2 - D[37]*Ts*u2_in))/Ts_e_tau2;
	P[36] = D[36];
	P[37] = D[37];
	P[38] = D[38] + Q[12];
	P[39] = D[39] + D[40]*Ts*e_b3 - D[43]*Ts*e_b3 - D[41]*Ts*e_b3*(bias3 - u3);
	P[40] = (e_tau_y*(D[40]*Ts + D[40]*e_tau_y + D[42]*Ts*u3 - D[42]*Ts*u3_in))/Ts_e_tau_y2;
	P[41] = D[41];
	P[42] = D[42];
	P[43] = D[43] + Q[13];

    
	/********* this is the update part of the equation ***********/

    float S[3] = {P[0] + s_a, P[1] + s_a, P[2] + s_a};

	X[0] = w1 + (2*P[0]*(gyro_x - w1))/S[0];
	X[1] = w2 + (2*P[1]*(gyro_y - w2))/S[1];
	X[2] = w3 + (2*P[2]*(gyro_z - w3))/S[2];
	X[3] = u1 + (2*P[3]*(gyro_x - w1))/S[0];
	X[4] = u2 + (2*P[5]*(gyro_y - w2))/S[1];
	X[5] = u3 + (2*P[7]*(gyro_z - w3))/S[2];
	X[6] = b1 + (2*P[9]*(gyro_x - w1))/S[0];
	X[7] = b2 + (2*P[12]*(gyro_y - w2))/S[1];
	X[8] = b3 + (2*P[15]*(gyro_z - w3))/S[2];
	X[9] = tau + (2*P[18]*(gyro_x - w1))/S[0] + (2*P[19]*(gyro_y - w2))/S[1];
	X[10] = tau_y + (2*P[25]*(gyro_z - w3))/S[2];
	X[11] = bias1 + (2*P[29]*(gyro_x - w1))/S[0];
	X[12] = bias2 + (2*P[34]*(gyro_y - w2))/S[1];
	X[13] = bias3 + (2*P[39]*(gyro_z - w3))/S[2];

	// update the duplicate cache
	for (uint32_t i = 0; i < RTSI_NUMP; i++)
        D[i] = P[i];
    
	// This is an approximation that removes some cross axis uncertainty but
	// substantially reduces the number of calculations
	// Covariance calculation
	P[0] = -D[0]*(D[0]/S[0] - 1);
	P[1] = -D[1]*(D[1]/S[1] - 1);
	P[2] = -D[2]*(D[2]/S[2] - 1);
	P[3] = -D[3]*(D[0]/S[0] - 1);
	P[4] = D[4] - D[3]*D[3]/S[0];
	P[5] = -D[5]*(D[1]/S[1] - 1);
	P[6] = D[6] - D[5]*D[5]/S[1];
	P[7] = -D[7]*(D[2]/S[2] - 1);
	P[8] = D[8] - D[7]*D[7]/S[2];
	P[9] = -D[9]*(D[0]/S[0] - 1);
	P[10] = D[10] - (D[3]*D[9])/S[0];
	P[11] = D[11] - D[9]*D[9]/S[0];
	P[12] = -D[12]*(D[1]/S[1] - 1);
	P[13] = D[13] - (D[5]*D[12])/S[1];
	P[14] = D[14] - D[12]*D[12]/S[1];
	P[15] = -D[15]*(D[2]/S[2] - 1);
	P[16] = D[16] - (D[7]*D[15])/S[2];
	P[17] = D[17] - D[15]*D[15]/S[2];
	P[18] = -D[18]*(D[0]/S[0] - 1);
	P[19] = -D[19]*(D[1]/S[1] - 1);
	P[20] = D[20] - (D[3]*D[18])/S[0];
	P[21] = D[21] - (D[5]*D[19])/S[1];
	P[22] = D[22] - (D[9]*D[18])/S[0];
	P[23] = D[23] - (D[12]*D[19])/S[1];
	P[24] = D[24] - D[18]*D[18]/S[0] - D[19]*D[19]/S[1];
	P[25] = -D[25]*(D[2]/S[2] - 1);
	P[26] = D[26] - (D[7]*D[25])/S[2];
	P[27] = D[27] - (D[15]*D[25])/S[2];
	P[28] = D[28] - D[25]*D[25]/S[2];
	P[29] = -D[29]*(D[0]/S[0] - 1);
	P[30] = D[30] - (D[3]*D[29])/S[0];
	P[31] = D[31] - (D[9]*D[29])/S[0];
	P[32] = D[32] - (D[18]*D[29])/S[0];
	P[33] = D[33] - D[29]*D[29]/S[0];
	P[34] = -D[34]*(D[1]/S[1] - 1);
	P[35] = D[35] - (D[5]*D[34])/S[1];
	P[36] = D[36] - (D[12]*D[34])/S[1];
	P[37] = D[37] - (D[19]*D[34])/S[1];
	P[38] = D[38] - D[34]*D[34]/S[1];
	P[39] = -D[39]*(D[2]/S[2] - 1);
	P[40] = D[40] - (D[7]*D[39])/S[2];
	P[41] = D[41] - (D[15]*D[39])/S[2];
	P[42] = D[42] - (D[25]*D[39])/S[2];
	P[43] = D[43] - D[39]*D[39]/S[2];


	// apply limits to some of the state variables
	if (X[9] > -1.5f)
	    X[9] = -1.5f;
	if (X[9] < -5.0f)
	    X[9] = -5.0f;
	if (X[10] > -1.5f)
	    X[10] = -1.5f;
	if (X[10] < -8.0f)
	    X[10] = -8.0f;
	if (X[11] > 0.5f)
	    X[11] = 0.5f;
	if (X[11] < -0.5f)
	    X[11] = -0.5f;
	if (X[12] > 0.5f)
	    X[12] = 0.5f;
	if (X[12] < -0.5f)
	    X[12] = -0.5f;
	if (X[13] > 0.5f)
	    X[13] = 0.5f;
	if (X[13] < -0.5f)
	    X[13] = -0.5f;
}

/**
 * Initialize the state variable and covariance matrix
 * for the system identification EKF
 */
void rtsi_init(float X[RTSI_NUMX], float P[RTSI_NUMP])
{
	const float q_init[RTSI_NUMX] = {
		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		0.05f, 0.05f, 0.005f,
		0.05f, 0.05f,
		0.05f, 0.05f, 0.05f
	};

	X[0] = X[1] = X[2] = 0.0f;    // assume no rotation
	X[3] = X[4] = X[5] = 0.0f;    // and no net torque
	X[6] = X[7]        = 10.0f;   // medium amount of strength
	X[8]               = 7.0f;    // yaw
	X[9] = X[10] = -4.0f;         // and 50 ms time scale
	X[11] = X[12] = X[13] = 0.0f; // zero bias

	// P initialization
	// Could zero this like: *P = *((float [AF_NUMP]){});
	P[0] = q_init[0];
	P[1] = q_init[1];
	P[2] = q_init[2];
	P[3] = 0.0f;
	P[4] = q_init[3];
	P[5] = 0.0f;
	P[6] = q_init[4];
	P[7] = 0.0f;
	P[8] = q_init[5];
	P[9] = 0.0f;
	P[10] = 0.0f;
	P[11] = q_init[6];
	P[12] = 0.0f;
	P[13] = 0.0f;
	P[14] = q_init[7];
	P[15] = 0.0f;
	P[16] = 0.0f;
	P[17] = q_init[8];
	P[18] = 0.0f;
	P[19] = 0.0f;
	P[20] = 0.0f;
	P[21] = 0.0f;
	P[22] = 0.0f;
	P[23] = 0.0f;
	P[24] = q_init[9];
	P[25] = 0.0f;
	P[26] = 0.0f;
	P[27] = 0.0f;
	P[28] = q_init[10];
	P[29] = 0.0f;
	P[30] = 0.0f;
	P[31] = 0.0f;
	P[32] = 0.0f;
	P[33] = q_init[11];
	P[34] = 0.0f;
	P[35] = 0.0f;
	P[36] = 0.0f;
	P[37] = 0.0f;
	P[38] = q_init[12];
	P[39] = 0.0f;
	P[40] = 0.0f;
	P[41] = 0.0f;
	P[42] = 0.0f;
	P[43] = q_init[13];
}

/**
 * @}
 * @}
 */
