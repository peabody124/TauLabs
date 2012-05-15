//
//  esc_log_buffer.c
//  OpenPilotOSX
//
//  Created by James Cotton on 12/2/11.
//  Copyright 2011 OpenPilot. All rights reserved.
//

#include "pios.h"

#define BACKBUFFER_SIZE 12048
uint8_t back_buf[BACKBUFFER_SIZE];
uint32_t back_buf_point = 0;

void esc_logger_init()
{
	back_buf_point = 0;
}

int32_t esc_logger_put(const uint8_t * data, int32_t length)
{
	if((back_buf_point + length) > sizeof(back_buf))
		return -1;
	
	memcpy(&back_buf[back_buf_point], data, length);
	back_buf_point += length;
	
	return 0;
}


uint8_t * esc_logger_getbuffer()
{
	return back_buf;
}

uint32_t esc_logger_getlength() {
	return back_buf_point;
}