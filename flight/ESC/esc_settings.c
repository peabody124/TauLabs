//
//  esc_settings.c
//  OpenPilotOSX
//
//  Created by James Cotton on 11/30/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#include <stdio.h>
#include "pios.h"
#include "esc_fsm.h"
#include <pios_board_info.h>

const struct esc_config default_config = {
	.max_dc_change = 0.1 * PIOS_ESC_MAX_DUTYCYCLE,
	.kp = 5, //0.0005 * PID_SCALE,
	.ki = 1, //0.0001 * PID_SCALE,
	.kff = 1.3e-4 * PID_SCALE,
	.kff2 = -0.05 * PID_SCALE,
	.ilim = 500,
	.min_dc = 0,
	.max_dc = 0.90 * PIOS_ESC_MAX_DUTYCYCLE,
	.initial_startup_speed = 100,
	.final_startup_speed = 400,
	.startup_current_target = 20,
	.commutation_phase = 23,
	.soft_current_limit = 3000, /* 10 mA per unit */
	.hard_current_limit = 2350,
	.commutation_offset_us = 0,
	.magic = ESC_CONFIG_MAGIC,
};

static int32_t settings_crc(const struct esc_config * settings) 
{
	uint8_t * p1 = (uint8_t *) settings;
	
	CRC_ResetDR();
    for (int32_t i = 0; i < sizeof(*settings); )
    {
        uint32_t value = 0;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 0;  else value |= 0x000000ff; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 8;  else value |= 0x0000ff00; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 16; else value |= 0x00ff0000; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 24; else value |= 0xff000000; i++;
		CRC_CalcCRC(value);
	}
	return CRC_GetCRC();
}

/**
 * Tries to load the settings from the EEPROM reserved for configuration
 * @param[out] settings populated with new settings if succeed, untouched if failure
 * @return < 0 for failure or 0 for success
 */
int32_t esc_settings_load(struct esc_config * settings) 
{
	if (pios_board_info_blob.magic != PIOS_BOARD_INFO_BLOB_MAGIC)
		return -1;
	
	int32_t eeprom_size = pios_board_info_blob.ee_size;
	
	if (eeprom_size < sizeof(struct esc_config))
		return -2;

	const struct esc_config * flash_settings = (const struct esc_config *) pios_board_info_blob.ee_base;
	if(flash_settings->magic != ESC_CONFIG_MAGIC)
		return -3;

	// Align and read the CRC
	uint32_t crc = settings_crc(flash_settings);
	uint32_t * crc_flash = (uint32_t *) (((uint32_t) flash_settings + sizeof(*settings) + 3) & 0xfffffffc);
	
	if(crc != *crc_flash)
		return -4;

	if(settings)
		memcpy((uint8_t *) settings, (uint8_t *) flash_settings, sizeof(*flash_settings));

	return 0;
}

int32_t esc_settings_defaults(struct esc_config * settings) 
{
	memcpy((uint8_t *) settings, (uint8_t *) &default_config, sizeof(default_config));
	
	return 0;
}

/**
 * Save the settings into EEPROM
 * @param[in] settings the settings to save
 * @returns 0 if success, -1 for failure
 */
int32_t esc_settings_save(struct esc_config * settings)
{   // save the settings into EEPROM
	
    FLASH_Status fs;
    uint32_t flash_addr;
    uint8_t *p1;
    uint32_t *p3;
	
	if (pios_board_info_blob.magic != PIOS_BOARD_INFO_BLOB_MAGIC)
		return -1;

    // address of settings in FLASH area
    flash_addr = pios_board_info_blob.ee_base;
	
	if(pios_board_info_blob.ee_size < sizeof(*settings))
	   return -2;
		
	uint32_t crc = settings_crc(settings);
	
    // Unlock the Flash Program Erase controller
    FLASH_Unlock();
	
	// TODO: Check the settings have change AND the area isn't already erased

	// Erase page
	fs = FLASH_ErasePage(pios_board_info_blob.ee_base);
	if (fs != FLASH_COMPLETE)
	{   // error
		FLASH_Lock();
		return -3;
	}	
	
    // save the settings into flash area (emulated EEPROM area)	
    p1 = (uint8_t *)  settings;
    p3 = (uint32_t *) flash_addr;
	
    // write 4 bytes at a time into program flash area (emulated EEPROM area)
    for (int32_t i = 0; i < sizeof(*settings); p3++)
    {
        uint32_t value = 0;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 0;  else value |= 0x000000ff; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 8;  else value |= 0x0000ff00; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 16; else value |= 0x00ff0000; i++;
        if (i < sizeof(*settings)) value |= (uint32_t)*p1++ << 24; else value |= 0xff000000; i++;
		
        fs = FLASH_ProgramWord((uint32_t)p3, value); // write a 32-bit value
        if (fs != FLASH_COMPLETE)
        {
            FLASH_Lock();			
            return -4;
		}
    }
	
	// Align and write the CRC
	p3 = (uint32_t *) ((flash_addr + sizeof(*settings) + 3) & 0xfffffffc);
	fs = FLASH_ProgramWord((uint32_t)p3, crc);
	if (fs != FLASH_COMPLETE)
		return -5;
	
    // Lock the Flash Program Erase controller
    FLASH_Lock();
	
	// Check that the values wolud be loadable
	return esc_settings_load(NULL) < 0 ? -6 : 0;
}
