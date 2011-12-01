//
//  esc_settings.h
//  OpenPilotOSX
//
//  Created by James Cotton on 11/30/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#ifndef OpenPilotOSX_esc_settings_h
#define OpenPilotOSX_esc_settings_h

#include "esc_fsm.h"

int32_t esc_settings_load(struct esc_config * settings);
int32_t esc_settings_defaults(struct esc_config * settings);
int32_t esc_settings_save(struct esc_config * settings);

#endif
