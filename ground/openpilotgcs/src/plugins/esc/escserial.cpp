#include "escserial.h"
#include "escsettings.h"
#include <QThread>
#include <QDebug>

/* The size of the data packets.  Unfortunately the compiler won't allow */
/* specified initializers so leaving comments in place to describe */
const int EscSerial::esc_command_data_size[ESC_COMMAND_LAST] = {
    /*[EscSerial::ESC_COMMAND_SET_CONFIG] = */ sizeof(EscSettings::DataFields) + 0, // No CRC yet
    /*[EscSerial::ESC_COMMAND_GET_CONFIG] = */ 0,
    /*[EscSerial::ESC_COMMAND_SAVE_CONFIG] = */ 0,
    /*[EscSerial::ESC_COMMAND_ENABLE_SERIAL_LOGGING] = */ 0,
    /*[EscSerial::ESC_COMMAND_DISABLE_SERIAL_LOGGING] = */ 0,
    /*[EscSerial::ESC_COMMAND_REBOOT_BL] = */ 0,
    /*[EscSerial::ESC_COMMAND_ENABLE_SERIAL_CONTROL] = */ 0,
    /*[EscSerial::ESC_COMMAND_DISABLE_SERIAL_CONTROL] = */ 0,
    /*[EscSerial::ESC_COMMAND_SET_SPEED] = */ 2,
    /*[EscSerial::ESC_COMMAND_WHOAMI] = */ 0,
    /*[EscSerial::ESC_COMMAND_ENABLE_ADC_LOGGING] = */ 0,
    /*[EscSerial::ESC_COMMAND_GET_ADC_LOG] = */ 0,
    /*[EscSerial::ESC_COMMAND_SET_PWM_FREQ] = */ 2,
    /*[EscSerial::ESC_COMMAND_GET_STATUS] = */ 0,
    /*[EscSerial::ESC_COMMAND_BOOTLOADER] = */ 2,
};

class SleepThread : public QThread
{
public:
    static void usleep(unsigned long usecs)
    {
        QThread::usleep(usecs);
    }
};

/**
  * Constructor
  * @param qio_in QIODevice to store
  */
EscSerial::EscSerial(QIODevice *qio_in) :
    qio(qio_in)
{
}

EscSerial::~EscSerial()
{
    qDebug() << "Closing serial port";
    qio->close();
    delete qio;
    qio = NULL;
}

EscStatus::DataFields EscSerial::getStatus()
{
    EscStatus::DataFields escStatusData;
    qint64 bytesRead;

    qio->readAll();
    writeCommand(ESC_COMMAND_GET_STATUS, NULL);
    bytesRead = qio->read((char *) &escStatusData, sizeof(escStatusData));

    return escStatusData;
}

EscSettings::DataFields EscSerial::getSettings()
{
    EscSettings::DataFields escSettingsData;
    qint64 bytesRead;

    qio->readAll();
    writeCommand(ESC_COMMAND_GET_CONFIG, NULL);
    bytesRead = qio->read((char *) &escSettingsData, sizeof(escSettingsData));

    return escSettingsData;
}

void EscSerial::setSettings(EscSettings::DataFields settings)
{
    writeCommand(ESC_COMMAND_SET_CONFIG, (char *) &settings);
}

void EscSerial::saveSettings()
{
    qDebug() << "Saving settings";
    writeCommand(ESC_COMMAND_SAVE_CONFIG, NULL);
}

void EscSerial::bootloader()
{
    char reboot[2] = {REBOOT_1,REBOOT_2};
    writeCommand(ESC_COMMAND_BOOTLOADER, reboot);
}

/**
  * Send command and associated payload
  */
void EscSerial::writeCommand(enum esc_serial_command command, const char *data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);

    // Prepare buffer of contents
    char * buffer = new char[2 + esc_command_data_size[command]];
    buffer[0] = ESC_SYNC_BYTE;
    buffer[1] = command;
    for (int i = 0; i < esc_command_data_size[command]; i++)
        buffer[i+2] = data[i];

    // Sending byte by byte because weird things happen otherwise
    for (int i = 0; i < 2 + esc_command_data_size[command]; i++) {
        qio->write(&buffer[i], 1);
        SleepThread::usleep(500);
    }
    delete buffer;
}

void EscSerial::readCommand(enum esc_serial_command command, char *data)
{
    Q_ASSERT(command < EscSerial::ESC_COMMAND_LAST);
    const char to_write[] = {ESC_SYNC_BYTE, command};
    qio->write(to_write,2);
    qio->read(data,esc_command_data_size[command]);

}
