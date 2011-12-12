/**
 ******************************************************************************
 * @addtogroup OpenPilotModules OpenPilot Modules
 * @{
 * @addtogroup Sensors
 * @brief Acquires sensor data 
 * Specifically updates the the @ref Gyros, @ref Accels, and @ref Magnetometer objects
 * @{
 *
 * @file       sensors.c
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @brief      Module to handle all comms to the AHRS on a periodic basis.
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

/**
 * Input objects: None, takes sensor data via pios
 * Output objects: @ref Gyros @ref Accels @ref Magnetometer
 *
 * The module executes in its own thread.
 *
 * UAVObjects are automatically generated by the UAVObjectGenerator from
 * the object definition XML file.
 *
 * Modules have no API, all communication to other modules is done through UAVObjects.
 * However modules may use the API exposed by shared libraries.
 * See the OpenPilot wiki for more details.
 * http://www.openpilot.org/OpenPilot_Application_Architecture
 *
 */

#include "pios.h"
#include "attitude.h"
#include "magnetometer.h"
#include "accels.h"
#include "gyros.h"
#include "gyrosbias.h"
#include "attitudeactual.h"
#include "attitudesettings.h"
#include "flightstatus.h"
#include "CoordinateConversions.h"

// Private constants
#define STACK_SIZE_BYTES 1540
#define TASK_PRIORITY (tskIDLE_PRIORITY+3)

#define F_PI 3.14159265358979323846f
#define PI_MOD(x) (fmod(x + F_PI, F_PI * 2) - F_PI)
// Private types

// Private variables
static xTaskHandle sensorsTaskHandle;

// Private functions
static void SensorsTask(void *parameters);
static void settingsUpdatedCb(UAVObjEvent * objEv);

static float gyroGain = 0.42;
static int16_t accelbias[3];
static float R[3][3];
static int8_t rotate = 0;
static bool zero_during_arming = false;

// These values are initialized by settings but can be updated by the attitude algorithm
static bool bias_correct_gyro = true;
static float gyro_bias[3] = {0,0,0};

/**
 * API for sensor fusion algorithms:
 * Configure(xQueueHandle gyro, xQueueHandle accel, xQueueHandle mag, xQueueHandle baro)
 *   Stores all the queues the algorithm will pull data from
 * FinalizeSensors() -- before saving the sensors modifies them based on internal state (gyro bias)
 * Update() -- queries queues and updates the attitude estiamte
 */


/**
 * Initialise the module.  Called before the start function
 * \returns 0 on success or -1 if initialisation failed
 */
int32_t SensorsInitialize(void)
{
	GyrosInitialize();
	GyrosBiasInitialize();
	AccelsInitialize();
	MagnetometerInitialize();
	AttitudeSettingsInitialize();
	
	for(uint8_t i = 0; i < 3; i++)
		for(uint8_t j = 0; j < 3; j++)
			R[i][j] = 0;
	
	AttitudeSettingsConnectCallback(&settingsUpdatedCb);
	
	return 0;
}

/**
 * Start the task.  Expects all objects to be initialized by this point.
 * \returns 0 on success or -1 if initialisation failed
 */
int32_t SensorsStart(void)
{
	// Start main task
	xTaskCreate(SensorsTask, (signed char *)"Sensors", STACK_SIZE_BYTES/4, NULL, TASK_PRIORITY, &sensorsTaskHandle);
	TaskMonitorAdd(TASKINFO_RUNNING_SENSORS, sensorsTaskHandle);
	PIOS_WDG_RegisterFlag(PIOS_WDG_SENSORS);

	return 0;
}

MODULE_INITCALL(SensorsInitialize, SensorsStart)

int32_t accel_test;
int32_t gyro_test;
int32_t mag_test;
//int32_t pressure_test;


/**
 * The sensor task.  This polls the gyros at 500 Hz and pumps that data to
 * stabilization and to the attitude loop
 */
static void SensorsTask(void *parameters)
{
	uint8_t init = 0;
	portTickType lastSysTime;
	uint32_t accel_samples;
	uint32_t gyro_samples;
	struct pios_bma180_data accel;
	struct pios_mpu6000_data gyro;
	int32_t accel_accum[3] = {0, 0, 0};
	int32_t gyro_accum[3] = {0,0,0};
	float scaling;

	AlarmsClear(SYSTEMALARMS_ALARM_SENSORS);

	accel_test = PIOS_BMA180_Test();
	gyro_test = PIOS_MPU6000_Test();
	mag_test = PIOS_HMC5883_Test();
	
	if(accel_test < 0 || gyro_test < 0 || mag_test < 0) {
		AlarmsSet(SYSTEMALARMS_ALARM_SENSORS, SYSTEMALARMS_ALARM_CRITICAL);
		while(1) {
			PIOS_WDG_UpdateFlag(PIOS_WDG_SENSORS);
			vTaskDelay(10);
		}
	}
	
	// Main task loop
	lastSysTime = xTaskGetTickCount();
	while (1) {
		// TODO: add timeouts to the sensor reads and set an error if the fail
	
		int32_t read_good;
		int32_t count;
		
		for (int i = 0; i < 3; i++) {
			accel_accum[i] = 0;
			gyro_accum[i] = 0;
		}
		accel_samples = 0;
		gyro_samples = 0;
		
		// Make sure we get one sample
		count = 0;
		while((read_good = PIOS_BMA180_ReadFifo(&accel)) != 0);
		while(read_good == 0) {	
			count++;
			
			accel_accum[0] += accel.x;
			accel_accum[1] += accel.y;
			accel_accum[2] += accel.z;
			
			read_good = PIOS_BMA180_ReadFifo(&accel);
		}
		accel_samples = count;
		
		float accels[3] = {(float) accel_accum[1] / accel_samples, (float) accel_accum[0] / accel_samples, -(float) accel_accum[2] / accel_samples};
		
		// Not the swaping of channel orders
		scaling = PIOS_BMA180_GetScale();
		AccelsData accelsData; // Skip get as we set all the fields
		accelsData.x = (accels[0] - accelbias[0]) * scaling;
		accelsData.y = (accels[1] - accelbias[1]) * scaling;
		accelsData.z = (accels[2] - accelbias[2]) * scaling;
		accelsData.temperature = 25.0f + ((float) accel.temperature - 2.0f) / 2.0f;
		AccelsSet(&accelsData);
		
		
		// Make sure we get one sample
		count = 0;
		while((read_good = PIOS_MPU6000_ReadFifo(&gyro)) != 0);
		while(read_good == 0) {
			count++;
			
			gyro_accum[0] += gyro.gyro_x;
			gyro_accum[1] += gyro.gyro_y;
			gyro_accum[2] += gyro.gyro_z;
			
			read_good = PIOS_MPU6000_ReadFifo(&gyro);
		}
		gyro_samples = count;
		
		float gyros[3] = {(float) gyro_accum[1] / gyro_samples, (float) gyro_accum[0] / gyro_samples, -(float) gyro_accum[2] / gyro_samples};
		
		scaling = PIOS_MPU6000_GetScale();
		GyrosData gyrosData; // Skip get as we set all the fields
		gyrosData.x = gyros[0] * scaling;
		gyrosData.y = gyros[1] * scaling;
		gyrosData.z = gyros[2] * scaling;
		gyrosData.temperature = 35.0f + ((float) gyro.temperature + 512.0f) / 340.0f;

		if (bias_correct_gyro) {
			// Apply bias correction to the gyros
			GyrosBiasData gyrosBias;
			GyrosBiasGet(&gyrosBias);
			gyrosData.x += gyrosBias.x;
			gyrosData.y += gyrosBias.y;
			gyrosData.z += gyrosBias.z;
		}
		
		GyrosSet(&gyrosData);
		
		// Because most crafts wont get enough information from gravity to zero yaw gyro, we try
		// and make it average zero (weakly)
		
		if (PIOS_HMC5883_NewDataAvailable()) {
			int16_t values[3];
			PIOS_HMC5883_ReadMag(values);
			MagnetometerData mag; // Skip get as we set all the fields
			mag.x = -values[0];
			mag.y = -values[1];
			mag.z = -values[2];
			MagnetometerSet(&mag);
		}
		
		PIOS_WDG_UpdateFlag(PIOS_WDG_SENSORS);
		vTaskDelayUntil(&lastSysTime, 2 / portTICK_RATE_MS);
	}
}

/**
 * Locally cache some variables from the AtttitudeSettings object
 */
static void settingsUpdatedCb(UAVObjEvent * objEv) {
	AttitudeSettingsData attitudeSettings;
	AttitudeSettingsGet(&attitudeSettings);

	gyroGain = attitudeSettings.GyroGain;

	zero_during_arming = attitudeSettings.ZeroDuringArming == ATTITUDESETTINGS_ZERODURINGARMING_TRUE;
	bias_correct_gyro = attitudeSettings.BiasCorrectGyro == ATTITUDESETTINGS_BIASCORRECTGYRO_TRUE;

	accelbias[0] = attitudeSettings.AccelBias[ATTITUDESETTINGS_ACCELBIAS_X];
	accelbias[1] = attitudeSettings.AccelBias[ATTITUDESETTINGS_ACCELBIAS_Y];
	accelbias[2] = attitudeSettings.AccelBias[ATTITUDESETTINGS_ACCELBIAS_Z];

	// Indicates not to expend cycles on rotation
	if(attitudeSettings.BoardRotation[0] == 0 && attitudeSettings.BoardRotation[1] == 0 &&
	   attitudeSettings.BoardRotation[2] == 0) {
		rotate = 0;

		// Shouldn't be used but to be safe
		float rotationQuat[4] = {1,0,0,0};
		Quaternion2R(rotationQuat, R);
	} else {
		float rotationQuat[4];
		const float rpy[3] = {attitudeSettings.BoardRotation[ATTITUDESETTINGS_BOARDROTATION_ROLL],
			attitudeSettings.BoardRotation[ATTITUDESETTINGS_BOARDROTATION_PITCH],
			attitudeSettings.BoardRotation[ATTITUDESETTINGS_BOARDROTATION_YAW]};
		RPY2Quaternion(rpy, rotationQuat);
		Quaternion2R(rotationQuat, R);
		rotate = 1;
	}
}
/**
  * @}
  * @}
  */
