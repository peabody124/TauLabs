#ifndef STM32_H
#define STM32_H

/*
  stm32flash - Open Source ST STM32 flash program for *nix
  Copyright (C) 2010 Geoffrey McRae <geoff@spacevs.com>

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <QIODevice>
#include <stdint.h>

typedef struct stm32     stm32_t;
typedef struct stm32_cmd	stm32_cmd_t;
typedef struct stm32_dev	stm32_dev_t;

struct stm32 {
        uint8_t			bl_version;
        uint8_t			version;
        uint8_t			option1, option2;
        uint16_t		pid;
        stm32_cmd_t		*cmd;
        const stm32_dev_t	*dev;
};

struct stm32_dev {
        uint16_t	id;
        char		*name;
        uint32_t	ram_start, ram_end;
        uint32_t	fl_start, fl_end;
        uint16_t	fl_pps; // pages per sector
        uint16_t	fl_ps;  // page size
        uint32_t	opt_start, opt_end;
        uint32_t	mem_start, mem_end;
};

class Stm32Bl {
public:
    Stm32Bl();
    ~Stm32Bl();

    stm32_t* stm32_init      (const char init);
    void stm32_close         ();
    char stm32_read_memory   (uint32_t address, uint8_t data[], unsigned int len);
    char stm32_write_memory  (uint32_t address, uint8_t data[], unsigned int len);
    char stm32_wunprot_memory();
    char stm32_erase_memory  (uint8_t spage, uint8_t pages);
    char stm32_go            (uint32_t address);
    char stm32_reset_device  ();

    int32_t uploadCode(QByteArray data);
    int32_t openDevice(QString serialPort);
    void print_device();
private:
    // Private variables
    QIODevice   *qio;
    stm32_t     *stm;
    double       uploaded;


    uint8_t stm32_gen_cs(const uint32_t v);
    void    stm32_send_byte(uint8_t byte);
    uint8_t stm32_read_byte();
    char    stm32_send_command(const uint8_t cmd);


    static const quint8 STM32_ACK =      0x79;
    static const quint8 STM32_NACK =     0x1F;
    static const quint8 STM32_CMD_INIT = 0x7F;
    static const quint8 STM32_CMD_GET =  0x00;	/* get the version and command supported */
    static const quint8 STM32_CMD_EE =   0x44;	/* extended erase */
};



#endif // STM32_H
