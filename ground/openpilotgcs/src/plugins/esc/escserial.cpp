#include "escserial.h"

//! The size of the data packets
const int esc_command_data_size[] = {
        [ESC_COMMAND_SET_CONFIG] = sizeof(EscSettingsData) + 0, // No CRC yet
        [ESC_COMMAND_GET_CONFIG] = 0,
        [ESC_COMMAND_SAVE_CONFIG] = 0,
        [ESC_COMMAND_ENABLE_SERIAL_LOGGING] = 0,
        [ESC_COMMAND_DISABLE_SERIAL_LOGGING] = 0,
        [ESC_COMMAND_REBOOT_BL] = 0,
        [ESC_COMMAND_DISABLE_SERIAL_CONTROL] = 0,
        [ESC_COMMAND_SET_SPEED] = 2,
        [ESC_COMMAND_SET_PWM_FREQ] = 2,
        [ESC_COMMAND_GET_STATUS] = 0,
};

/**
  * Constructor
  * @param qio_in QIODevice to store
  */
EscSerial::EscSerial(QIODevice qio_in) :
    qio(qio_in)
{
}

void EscSerial::getStatus()
{

}

void EscSerial::setSettings()
{

}

void EscSerial::writeCommand(enum esc_serial_command command, const char * data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);
    const char to_write[] = {ESC_SYNC_BYTE, command};
    qio.writeData(to_write,2);
    qio.write(data,esc_command_data_size[command]);
}

void EscSerial::readCommand(enum esc_serial_command command, const char * data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);
}
