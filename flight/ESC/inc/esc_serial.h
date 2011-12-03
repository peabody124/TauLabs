//
//  esc_serial.h
//  OpenPilotOSX
//
//  Created by James Cotton on 12/1/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#ifndef OpenPilotOSX_esc_serial_h
#define OpenPilotOSX_esc_serial_h

int32_t esc_serial_init();
int32_t esc_serial_parse(int32_t c);
int32_t esc_serial_process();

#endif
