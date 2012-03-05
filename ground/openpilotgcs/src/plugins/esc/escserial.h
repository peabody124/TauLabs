#ifndef ESCSERIAL_H
#define ESCSERIAL_H

#include "escstatus.h"
#include "escsettings.h"

class EscSerial
{
public:
    EscSerial(QIODevice *qio);
    ~EscSerial();

    void saveSettings();
    EscStatus::DataFields getStatus();
    EscSettings::DataFields getSettings();
    void setSettings(EscSettings::DataFields);

    void bootloader();

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
            ESC_COMMAND_BOOTLOADER = 14,
            ESC_COMMAND_LAST = 15
    };

    static const int esc_command_data_size[ESC_COMMAND_LAST];
    static const quint8 REBOOT_1 = 0x73;
    static const quint8 REBOOT_2 = 0x37;
private:

    static const quint8 ESC_SYNC_BYTE = 0x85;

    void writeCommand(enum esc_serial_command command, const char *data);
    void readCommand(enum esc_serial_command command, char *data);

    EscSettings::DataFields escSettings;
    EscStatus::DataFields escStatus;
    QIODevice *qio;
};

#endif // ESCSERIAL_H
