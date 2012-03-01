#include "escserial.h"
#include "escsettings.h"

/* The size of the data packets.  Unfortunately the compiler won't allow */
/* specified initializers so leaving comments in place to describe */
const int EscSerial::esc_command_data_size[ESC_COMMAND_LAST] = {
        /*[EscSerial::ESC_COMMAND_SET_CONFIG] = */ sizeof(EscSettings::DataFields) + 0, // No CRC yet
        /*[EscSerial::ESC_COMMAND_GET_CONFIG] = */ 0,
        /*[EscSerial::ESC_COMMAND_SAVE_CONFIG] = */ 0,
        /*[EscSerial::ESC_COMMAND_ENABLE_SERIAL_LOGGING] = */ 0,
        /*[EscSerial::ESC_COMMAND_DISABLE_SERIAL_LOGGING] = */ 0,
        /*[EscSerial::ESC_COMMAND_REBOOT_BL] = */ 0,
        /*[EscSerial::ESC_COMMAND_DISABLE_SERIAL_CONTROL] = */ 0,
        /*[EscSerial::ESC_COMMAND_SET_SPEED] = */ 2,
        /*[EscSerial::ESC_COMMAND_SET_PWM_FREQ] = */ 2,
        /*[EscSerial::ESC_COMMAND_GET_STATUS] = */ 0,
};

/**
  * Constructor
  * @param qio_in QIODevice to store
  */
EscSerial::EscSerial(QIODevice *qio_in) :
    qio(qio_in)
{
}

void EscSerial::getStatus()
{

}

void EscSerial::setSettings()
{

}

void EscSerial::writeCommand(enum esc_serial_command command, const char *data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);
    const char to_write[] = {ESC_SYNC_BYTE, command};
    qio->write(to_write,2);
    qio->write(data,esc_command_data_size[command]);
}

void EscSerial::readCommand(enum esc_serial_command command, char *data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);
    const char to_write[] = {ESC_SYNC_BYTE, command};
    qio->write(to_write,2);
    qio->read(data,esc_command_data_size[command]);

}
