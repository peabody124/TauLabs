//
//  esc_settings.h
//  OpenPilotOSX
//
//  Created by James Cotton on 11/30/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#ifndef OpenPilotOSX_esc_settings_h
#define OpenPilotOSX_esc_settings_h

typedef void *xQueueHandle;
#include "uavobjectmanager.h"

#include "escsettings.h"
#include "escstatus.h"

int32_t esc_settings_load(EscSettingsData * settings);
int32_t esc_settings_defaults(EscSettingsData * settings);
int32_t esc_settings_save(EscSettingsData * settings);

#endif
