//
//  esc_logger.h
//  OpenPilotOSX
//
//  Created by James Cotton on 12/2/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#ifndef OpenPilotOSX_esc_logger_h
#define OpenPilotOSX_esc_logger_h

void esc_logger_init();
int32_t esc_logger_put(const uint8_t * data, int32_t length);
uint8_t * esc_logger_getbuffer();
uint32_t esc_logger_getlength();

#endif
