/**
 ******************************************************************************
 * @addtogroup TauLabsModules Tau Labs Modules
 * @{
 * @addtogroup Rate Torque Linear Quadratic Gaussian controller
 * @{
 *
 * @file       rate_torque_lqr.h
 * @author     Tau Labs, http://taulabs.org, Copyright (C) 2016
 * @brief      LQR controller based on rate and torque estimate
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

#ifndef RATE_TORQUE_LQR_H
#define RATE_TORQUE_LQR_H

#include "stdint.h"
 
float rtlqr_calculate_axis(const float *rtkf_X, float rate_desired, uint32_t axis);

#endif /* RATE_TORQUE_LQR_H */
 /**
 * @}
 * @}
 */
 