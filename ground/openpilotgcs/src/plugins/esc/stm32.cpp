
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

#include <QtEndian>
#include <QDebug>
#include <qextserialport.h>

#include "stm32.h"

struct stm32_cmd {
    uint8_t get;
    uint8_t gvr;
    uint8_t gid;
    uint8_t rm;
    uint8_t go;
    uint8_t wm;
    uint8_t er; /* this may be extended erase */
    uint8_t wp;
    uint8_t uw;
    uint8_t rp;
    uint8_t ur;
};

/* device table */
const stm32_dev_t devices[] = {
    {0x412, "Low-density"       , 0x20000200, 0x20002800, 0x08000000, 0x08008000,  4, 1024, 0x1FFFF800, 0x1FFFF80F, 0x1FFFF000, 0x1FFFF800},
    {0x410, "Medium-density"    , 0x20000200, 0x20005000, 0x08000000, 0x08020000,  4, 1024, 0x1FFFF800, 0x1FFFF80F, 0x1FFFF000, 0x1FFFF800},
    {0x411, "STM32F2xx"         , 0x20002000, 0x20020000, 0x08000000, 0x08100000,  4, 16384, 0x1FFFC000, 0x1FFFC00F, 0x1FFF0000, 0x1FFF77DF},
    {0x413, "STM32F4xx"         , 0x20002000, 0x20020000, 0x08000000, 0x08100000,  4, 16384, 0x1FFFC000, 0x1FFFC00F, 0x1FFF0000, 0x1FFF77DF},
    {0x414, "High-density"      , 0x20000200, 0x20010000, 0x08000000, 0x08080000,  2, 2048, 0x1FFFF800, 0x1FFFF80F, 0x1FFFF000, 0x1FFFF800},
    {0x416, "Medium-density ULP", 0x20000800, 0x20004000, 0x08000000, 0x08020000, 16,  256, 0x1FF80000, 0x1FF8000F, 0x1FF00000, 0x1FF01000},
    {0x418, "Connectivity line" , 0x20001000, 0x20010000, 0x08000000, 0x08040000,  2, 2048, 0x1FFFF800, 0x1FFFF80F, 0x1FFFB000, 0x1FFFF800},
    {0x420, "Medium-density VL" , 0x20000200, 0x20002000, 0x08000000, 0x08020000,  4, 1024, 0x1FFFF800, 0x1FFFF80F, 0x1FFFF000, 0x1FFFF800},
    {0x428, "High-density VL"   , 0x20000200, 0x20008000, 0x08000000, 0x08080000,  2, 2048, 0x1FFFF800, 0x1FFFF80F, 0x1FFFF000, 0x1FFFF800},
    {0x430, "XL-density"        , 0x20000800, 0x20018000, 0x08000000, 0x08100000,  2, 2048, 0x1FFFF800, 0x1FFFF80F, 0x1FFFE000, 0x1FFFF800},
    {0x0}
};

/**
  * This is a program which is uploaded and resets the devices
  */
const unsigned int stmreset_length = 428;
const unsigned char stmreset_binary[] = {
0x00,0x62,0x00,0x20,0x33,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x00,0x00,0x00,0x00,
0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x00,0x00,0x00,0x00,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,0x30,0x03,0x00,0x20,
0x70,0x47,0xdf,0xf8,0x3c,0x10,0x0d,0x4a,0x0e,0x4b,0x5b,0x1a,0x05,0xd0,0x12,0xf8,
0x01,0x4b,0x01,0xf8,0x01,0x4b,0x01,0x3b,0xf9,0xdc,0x0b,0x49,0x0b,0x4b,0x5b,0x1a,
0x05,0xd0,0x4f,0xf0,0x00,0x02,0x01,0xf8,0x01,0x2b,0x01,0x3b,0xfb,0xdc,0x4f,0xf0,
0x00,0x00,0x4f,0xf0,0x00,0x01,0x00,0xf0,0x0b,0xf8,0xfe,0xe7,0xac,0x03,0x00,0x20,
0x00,0x22,0x00,0x20,0x00,0x22,0x00,0x20,0x00,0x22,0x00,0x20,0x00,0x22,0x00,0x20,
0x80,0xb4,0x00,0xaf,0x4e,0xf6,0x00,0x52,0xce,0xf2,0x00,0x02,0x4e,0xf6,0x00,0x53,
0xce,0xf2,0x00,0x03,0xdb,0x68,0x03,0xf4,0xe0,0x61,0x40,0xf2,0x04,0x03,0xc0,0xf2,
0xfa,0x53,0x41,0xea,0x03,0x03,0xd3,0x60,0xfe,0xe7,0x00,0xbf};

/* detect CPU endian */
char cpu_le() {
        const uint32_t cpu_le_test = 0x12345678;
        return ((unsigned char*)&cpu_le_test)[0] == 0x78;
}

uint32_t be_u32(const uint32_t v) {
        if (cpu_le())
                return  ((v & 0xFF000000) >> 24) |
                        ((v & 0x00FF0000) >>  8) |
                        ((v & 0x0000FF00) <<  8) |
                        ((v & 0x000000FF) << 24);
        return v;
}

uint8_t Stm32Bl::stm32_gen_cs(const uint32_t v) {
    return  ((v & 0xFF000000) >> 24) ^
            ((v & 0x00FF0000) >> 16) ^
            ((v & 0x0000FF00) >>  8) ^
            ((v & 0x000000FF) >>  0);
}

Stm32Bl::Stm32Bl() : qio(NULL), stm(NULL), uploaded(0)
{
}



Stm32Bl::~Stm32Bl()
{
    stm32_close();
}

/**
  * @brief Open the serial port and initialize bootloader
  */
int Stm32Bl::openDevice(QString dev)
{

    PortSettings settings;
    settings.BaudRate = BAUD57600;
    settings.DataBits = DATA_8;
    settings.Parity = PAR_EVEN;
    settings.StopBits = STOP_1;
    settings.FlowControl = FLOW_OFF;
    settings.Timeout_Millisec = 1000;

    QextSerialPort *serial_dev = new QextSerialPort(dev, settings);
    if (!serial_dev) {
        qDebug() << "Failed to initial serial port";
        return -1;
    }

    if (!serial_dev->open(QIODevice::ReadWrite))
    {
        qDebug() << "Failed to open serial port";
        delete serial_dev;
        return -1;
    }


    qio = serial_dev;
    stm = stm32_init(1);

    if(stm==NULL) {
        qDebug() << "Failed to initialize bootloader";
        serial_dev->close();
        return -1;
    }

    return 0;
}

void Stm32Bl::stm32_send_byte(uint8_t byte) {
    Q_UNUSED(stm);

    qint64 bytes = qio->write((const char *) &byte,1);
    if (bytes != 1)
        qDebug() << "Failed to write";
    //qDebug() << QString("Sent %1").arg(byte);
}

uint8_t Stm32Bl::stm32_read_byte() {
    uint8_t byte;
    qint64 bytes;
    bytes = qio->read((char *) &byte,1);
    if (bytes != 1)
        qDebug() << "Failed to read";
    //else
    //    qDebug() << QString("Read %1").arg(byte);
    return byte;
}

char Stm32Bl::stm32_send_command(const uint8_t cmd) {
    stm32_send_byte(cmd);
    stm32_send_byte(cmd ^ 0xFF);
    if (stm32_read_byte() != STM32_ACK) {
        qDebug() << QString("Error sending command 0x%02x to device\n").arg(cmd);
        return 0;
    }
    return 1;
}

stm32_t* Stm32Bl::stm32_init(const char init) {
    uint8_t len;

    qDebug() << "Initializing bootloader";
    stm      = new stm32_t;
    stm->cmd = new stm32_cmd_t;

    if (init) {
        stm32_send_byte(STM32_CMD_INIT);
        if (stm32_read_byte() != STM32_ACK) {
            stm32_close();
            qDebug() << "Failed to get init ACK from device";
            return NULL;
        }
    }

    /* get the bootloader information */
    if (!stm32_send_command(STM32_CMD_GET)) return 0;
    len              = stm32_read_byte() + 1;
    stm->bl_version  = stm32_read_byte(); --len;
    stm->cmd->get    = stm32_read_byte(); --len;
    stm->cmd->gvr    = stm32_read_byte(); --len;
    stm->cmd->gid    = stm32_read_byte(); --len;
    stm->cmd->rm     = stm32_read_byte(); --len;
    stm->cmd->go     = stm32_read_byte(); --len;
    stm->cmd->wm     = stm32_read_byte(); --len;
    stm->cmd->er     = stm32_read_byte(); --len;
    stm->cmd->wp     = stm32_read_byte(); --len;
    stm->cmd->uw     = stm32_read_byte(); --len;
    stm->cmd->rp     = stm32_read_byte(); --len;
    stm->cmd->ur     = stm32_read_byte(); --len;
    if (len > 0) {
        fprintf(stderr, "Seems this bootloader returns more then we understand in the GET command, we will skip the unknown bytes\n");
        while(len-- > 0) stm32_read_byte();
    }
    if (stm32_read_byte() != STM32_ACK) {
        stm32_close();
        return NULL;
    }

    /* get the version and read protection status  */
    if (!stm32_send_command(stm->cmd->gvr)) {
        stm32_close();
        return NULL;
    }

    stm->version = stm32_read_byte();
    stm->option1 = stm32_read_byte();
    stm->option2 = stm32_read_byte();
    if (stm32_read_byte() != STM32_ACK) {
        stm32_close();
        return NULL;
    }

    /* get the device ID */
    if (!stm32_send_command(stm->cmd->gid)) {
        stm32_close();
        return NULL;
    }
    len = stm32_read_byte() + 1;
    if (len != 2) {
        stm32_close();
        qDebug() << "More then two bytes sent in the PID, unknown/unsupported device";
        return NULL;
    }
    stm->pid = (stm32_read_byte() << 8) | stm32_read_byte();
    if (stm32_read_byte() != STM32_ACK) {
        stm32_close();
        return NULL;
    }

    stm->dev = devices;
    while(stm->dev->id != 0x00 && stm->dev->id != stm->pid)
        ++stm->dev;

    if (!stm->dev->id) {
        qDebug() << QString("Unknown/unsupported device (Device ID: 0x%03x)").arg(stm->pid);
        stm32_close();
        return NULL;
    }

    return stm;
}

void Stm32Bl::stm32_close()
{
    if (stm) free(stm->cmd);
    free(stm);
    stm = NULL;
}

void Stm32Bl::print_device()
{
    printf("Version      : 0x%02x\n", stm->bl_version);
    printf("Option 1     : 0x%02x\n", stm->option1);
    printf("Option 2     : 0x%02x\n", stm->option2);
    printf("Device ID    : 0x%04x (%s)\n", stm->pid, stm->dev->name);
    printf("RAM          : %dKiB  (%db reserved by bootloader)\n", (stm->dev->ram_end - 0x20000000) / 1024, stm->dev->ram_start - 0x20000000);
    printf("Flash        : %dKiB (sector size: %dx%d)\n", (stm->dev->fl_end - stm->dev->fl_start ) / 1024, stm->dev->fl_pps, stm->dev->fl_ps);
    printf("Option RAM   : %db\n", stm->dev->opt_end - stm->dev->opt_start);
    printf("System RAM   : %dKiB\n", (stm->dev->mem_end - stm->dev->mem_start) / 1024);
}

char Stm32Bl::stm32_read_memory(uint32_t address, uint8_t data[], unsigned int len) {
    uint8_t cs;
    unsigned int i;
    qint64 bytes;
    Q_ASSERT(len > 0 && len < 257);

    /* must be 32bit aligned */
    Q_ASSERT(address % 4 == 0);

    address = be_u32(address);
    cs      = stm32_gen_cs(address);

    if (!stm32_send_command(stm->cmd->rm)) return 0;
    bytes = qio->write((const char*) &address,4);
    Q_ASSERT(bytes == 4);
    stm32_send_byte(cs);
    if (stm32_read_byte() != STM32_ACK) return 0;

    i = len - 1;
    stm32_send_byte(i);
    stm32_send_byte(i ^ 0xFF);
    if (stm32_read_byte() != STM32_ACK) return 0;

    bytes = qio->read((char *) data, len);
    if (bytes != len) {
        qDebug() << "Failed to read memory";
        return 0;
    }

    return 1;
}

char Stm32Bl::stm32_write_memory(uint32_t address, uint8_t data[], unsigned int len) {
    uint8_t cs;
    unsigned int i;
    int c, extra;
    qint64 bytes;
    Q_ASSERT(len > 0 && len < 257);

    /* must be 32bit aligned */
    Q_ASSERT(address % 4 == 0);

    address = be_u32(address);
    cs      = stm32_gen_cs(address);

    /* send the address and checksum */
    if (!stm32_send_command(stm->cmd->wm)) return 0;
    qio->write((const char* ) &address, 4);
    qDebug() << QString().sprintf("Writing to 0x%08x",address);
    stm32_send_byte(cs);
    if (stm32_read_byte() != STM32_ACK) return 0;

    /* setup the cs and send the length */
    extra = len % 4;
    cs = len - 1 + extra;
    stm32_send_byte(cs);

    /* write the data and build the checksum */
    for(i = 0; i < len; ++i)
        cs ^= data[i];

    bytes = qio->write((const char *) data, len);
    Q_ASSERT(bytes == len);

    /* write the alignment padding */
    for(c = 0; c < extra; ++c) {
        stm32_send_byte(0xFF);
        cs ^= 0xFF;
    }

    /* send the checksum */
    stm32_send_byte(cs);
    return stm32_read_byte() == STM32_ACK;
}

char Stm32Bl::stm32_wunprot_memory() {
    if (!stm32_send_command(stm->cmd->uw)) return 0;
    if (!stm32_send_command(0x8C        )) return 0;
    return 1;
}

char Stm32Bl::stm32_erase_memory(uint8_t spage, uint8_t pages) {
    if (!stm32_send_command(stm->cmd->er)) {
        qDebug() << "Can't initiate chip erase!\n";
        return 0;
    }

    /* The erase command reported by the bootloader is either 0x43 or 0x44 */
    /* 0x44 is Extended Erase, a 2 byte based protocol and needs to be handled differently. */
    if (stm->cmd->er == STM32_CMD_EE) {
        /* Not all chips using Extended Erase support mass erase */
        /* Currently known as not supporting mass erase is the Ultra Low Power STM32L15xx range */
        /* So if someone has not overridden the default, but uses one of these chips, take it out of */
        /* mass erase mode, so it will be done page by page. This maximum might not be correct either! */
        if (stm->pid == 0x416 && pages == 0xFF) pages = 0xF8; /* works for the STM32L152RB with 128Kb flash */

        if (pages == 0xFF) {
            stm32_send_byte(0xFF);
            stm32_send_byte(0xFF); // 0xFFFF the magic number for mass erase
            stm32_send_byte(0x00); // 0x00 the XOR of those two bytes as a checksum
            if (stm32_read_byte() != STM32_ACK) {
                qDebug() << "Mass erase failed. Try specifying the number of pages to be erased.\n";
                return 0;
            }
            return 1;
        }

        uint16_t pg_num;
        uint8_t pg_byte;
        uint8_t cs = 0;

        stm32_send_byte(pages >> 8); // Number of pages to be erased, two bytes, MSB first
        stm32_send_byte(pages & 0xFF);

        for (pg_num = 0; pg_num <= pages; pg_num++) {
            pg_byte = pg_num >> 8;
            cs ^= pg_byte;
            stm32_send_byte(pg_byte);
            pg_byte = pg_num & 0xFF;
            cs ^= pg_byte;
            stm32_send_byte(pg_byte);
        }
        stm32_send_byte(0x00);  // Ought to need to hand over a valid checksum here...but 0 seems to work!

        if (stm32_read_byte() != STM32_ACK) {
            qDebug() << "Page-by-page erase failed. Check the maximum pages your device supports.\n";
            return 0;
        }

        return 1;
    }

    /* And now the regular erase (0x43) for all other chips */
    if (pages == 0xFF) {
        return stm32_send_command(0xFF);
    } else {
        uint8_t cs = 0;
        uint8_t pg_num;
        stm32_send_byte(pages-1);
        cs ^= (pages-1);
        for (pg_num = spage; pg_num < (pages + spage); pg_num++) {
            stm32_send_byte(pg_num);
            cs ^= pg_num;
        }
        stm32_send_byte(cs);
        return stm32_read_byte() == STM32_ACK;
    }
}

/**
  * Jump to specified address or to flash start if zero
  * @param[in] address address to jump to
  */
char Stm32Bl::stm32_go(uint32_t address) {
    uint8_t cs;

    if(address == 0)
        address = stm->dev->fl_start;

    address = qToBigEndian(address);
    cs      = stm32_gen_cs(address);

    if (!stm32_send_command(stm->cmd->go)) return 0;
    qio->write((char *) &address, 4);
    qio->write((char *) &cs, 1);

    return stm32_read_byte() == STM32_ACK;
}

char Stm32Bl::stm32_reset_device() {
    /*
                since the bootloader does not have a reset command, we
                upload the stmreset program into ram and run it, which
                resets the device for us
        */

    uint32_t length		= stmreset_length;
    unsigned char* pos	= (unsigned char *) stmreset_binary;
    uint32_t address	= stm->dev->ram_start;
    while(length > 0) {
        uint32_t w = length > 256 ? 256 : length;
        if (!stm32_write_memory(address, pos, w))
            return 0;

        address	+= w;
        pos	+= w;
        length	-= w;
    }

    return stm32_go(stm->dev->ram_start);
}

/**
  * Upload new code to the ESC
  * @param[in] data code to upload
  * @return -1 for fail, 0 for success
  */
int32_t Stm32Bl::uploadCode(QByteArray data)
{
    off_t   offset = 0;
    ssize_t r;
    int spage = 0;
    int npages = 0xff;
    int len;
    quint32 addr;
    unsigned int size = data.length();

    if (size > stm->dev->fl_end - stm->dev->fl_start) {
        qDebug() << "File provided larger then available flash space.";
        return -1;
    }

    qDebug() << "Erasing memory";
    if (!stm32_erase_memory(spage, npages)) {
        qDebug() << "Failed to erase memory";
        return -1;
    }
    qDebug() << "Erasing memory done";

    uint8_t *buffer = (uint8_t *) data.data();

    int idx = 0;

    addr = stm->dev->fl_start + (spage * stm->dev->fl_ps);

    while(addr < stm->dev->fl_end && offset < size) {
        quint32 left   = stm->dev->fl_end - addr;
        len             = 256 > left ? left : 256;
        len             = len > size - offset ? size - offset : len;

        qDebug() << "Write memory";
        if (!stm32_write_memory(addr, buffer, len)) {
            qDebug() << QString().sprintf("Failed to write memory at address 0x%08x",addr);
            return -1;
        }
        qDebug() << "Write memory done";

        addr    += len;
        offset  += len;
        buffer  += len;

        uploaded = offset / (stm->dev->fl_end - stm->dev->fl_start);

        qDebug() << QString().sprintf("Wrote address 0x%08x, completed %f",addr, uploaded);
    }
}
