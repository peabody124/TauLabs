#ifndef ESCSERIAL_H
#define ESCSERIAL_H

#include "escstatus.h"
#include "escsettings.h"

class EscSerial
{
public:
    EscSerial(QIODevice qio);

    void getStatus();
    void setSettings();

private:
    enum esc_serial_command {
            ESC_COMMAND_SET_CONFIG = 0,
            ESC_COMMAND_GET_CONFIG = 1,
            ESC_COMMAND_SAVE_CONFIG = 2,
            ESC_COMMAND_ENABLE_SERIAL_LOGGING = 3,
            ESC_COMMAND_DISABLE_SERIAL_LOGGING = 4,
            ESC_COMMAND_REBOOT_BL = 5,
            ESC_COMMAND_ENABLE_SERIAL_CONTROL = 6,
            ESC_COMMAND_DISABLE_SERIAL_CONTROL = 7,
            ESC_COMMAND_SET_SPEED = 8,
            ESC_COMMAND_WHOAMI = 9,
            ESC_COMMAND_ENABLE_ADC_LOGGING = 10,
            ESC_COMMAND_GET_ADC_LOG = 11,
            ESC_COMMAND_SET_PWM_FREQ = 12,
            ESC_COMMAND_GET_STATUS = 13,
            ESC_COMMAND_LAST = 14
    };

    const quint8 ESC_SYNC_BYTE = 0x85;

    void writeCommand(enum esc_serial_command command, quint8 * data);
    void readCommand(enum esc_serial_command command, quint8 * data);

    EscSettings::DataFields escSettings;
    EscStatus::DataFields escStatus;
    QIODevice qio;
};

#endif // ESCSERIAL_H
