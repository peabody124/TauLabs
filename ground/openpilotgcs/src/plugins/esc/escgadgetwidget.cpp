/**
 ******************************************************************************
 *
 * @file       escgadgetwidget.cpp
 * @author     The OpenPilot Team, http://www.openpilot.org Copyright (C) 2010.
 * @addtogroup GCSPlugins GCS Plugins
 * @{
 * @{
 *****************************************************************************/
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <QDebug>
#include <QtOpenGL/QGLWidget>
#include <QSettings>
#include <QFileInfo>
#include <QDir>
#include <QFileDialog>
#include "escsettings.h"
#include <stdlib.h>

#include "escgadgetwidget.h"

#define NO_PORT         0
#define SERIAL_PORT     1
#define USB_PORT        2

// constructor
EscGadgetWidget::EscGadgetWidget(QWidget *parent) :
    QWidget(parent),
    m_widget(NULL),
    escSerial(NULL),
    connectedStatus(false)
{
    m_widget = new Ui_EscWidget();
    m_widget->setupUi(this);

    // Listen to telemetry connection events
    ExtensionSystem::PluginManager *pluginManager = ExtensionSystem::PluginManager::instance();
    if (pluginManager)
    {
        TelemetryManager *telemetryManager = pluginManager->getObject<TelemetryManager>();
        if (telemetryManager)
        {
            connect(telemetryManager, SIGNAL(myStart()), this, SLOT(onTelemetryStart()));
            connect(telemetryManager, SIGNAL(myStop()), this, SLOT(onTelemetryStop()));
            connect(telemetryManager, SIGNAL(connected()), this, SLOT(onTelemetryConnect()));
            connect(telemetryManager, SIGNAL(disconnected()), this, SLOT(onTelemetryDisconnect()));
        }
    }

    getPorts();

    setEnabled(true);

    m_widget->saveSettingsToFlash->setEnabled(false);
    m_widget->saveSettingsToRAM->setEnabled(false);

    connect(m_widget->connectButton, SIGNAL(clicked()), this, SLOT(connectDisconnect()));
    connect(m_widget->saveSettingsToRAM, SIGNAL(clicked()), this, SLOT(applyConfiguration()));
    connect(m_widget->saveSettingsToFlash, SIGNAL(clicked()), this, SLOT(saveConfiguration()));
}

// destructor .. this never gets called :(
EscGadgetWidget::~EscGadgetWidget()
{
}

void EscGadgetWidget::resizeEvent(QResizeEvent *event)
{
    Q_UNUSED(event);
}

void EscGadgetWidget::onComboBoxPorts_currentIndexChanged(int index)
{
    Q_UNUSED(index);
}

bool sortSerialPorts(const QextPortInfo &s1, const QextPortInfo &s2)
{
    return (s1.portName < s2.portName);
}

void EscGadgetWidget::getPorts()
{
    disconnect(m_widget->comboBox_Ports, 0, 0, 0);

    m_widget->comboBox_Ports->clear();

    // Add item for connecting via FC
    m_widget->comboBox_Ports->addItem("FlightController");

    // Populate the telemetry combo box with serial ports
    QList<QextPortInfo> serial_ports = QextSerialEnumerator::getPorts();
    qSort(serial_ports.begin(), serial_ports.end(), sortSerialPorts);
    QStringList list;
    foreach (QextPortInfo port, serial_ports)
        m_widget->comboBox_Ports->addItem("com: " + port.friendName, SERIAL_PORT);

    connect(m_widget->comboBox_Ports, SIGNAL(currentIndexChanged(int)), this, SLOT(onComboBoxPorts_currentIndexChanged(int)));

    onComboBoxPorts_currentIndexChanged(m_widget->comboBox_Ports->currentIndex());

}

void EscGadgetWidget::onTelemetryStart()
{

}

void EscGadgetWidget::onTelemetryStop()
{

}

void EscGadgetWidget::onTelemetryConnect()
{
}

void EscGadgetWidget::onTelemetryDisconnect()
{
}

void EscGadgetWidget::disableTelemetry()
{
    // Suspend telemetry & polling
}

void EscGadgetWidget::enableTelemetry()
{	// Restart the polling thread
}

void EscGadgetWidget::getSettings()
{
    const float PID_SCALE = 32178;

    EscSettings::DataFields escSettings = EscSettings::GetInstance(getObjectManager())->getData();

    // Check the UI elements
    escSettings.Mode = m_widget->closedLoopCheck->checkState() ? EscSettings::MODE_CLOSED : EscSettings::MODE_OPEN;
    escSettings.RisingKp = m_widget->spinRisingKp->value();
    escSettings.FallingKp = m_widget->spinFallingKp->value();
    escSettings.CommutationPhase = m_widget->spinPhase->value();
    escSettings.Kv = m_widget->spinKv->value();

    if(escSettings.RisingKp > PID_SCALE)
        escSettings.RisingKp = PID_SCALE;
    if(escSettings.FallingKp > PID_SCALE)
        escSettings.FallingKp = PID_SCALE;
}

void EscGadgetWidget::sendSettings()
{

}

QString EscGadgetWidget::getSerialPortDevice(const QString &friendName)
{
        QList<QextPortInfo> ports = QextSerialEnumerator::getPorts();

        foreach (QextPortInfo port, ports)
        {
                #ifdef Q_OS_WIN
                        if (port.friendName == friendName)
                                return port.portName;
                #else
                        if (port.friendName == friendName)
                                return port.physName;
                #endif
        }

        return "";
}

/**
  * @brief Disconnect the selected port
  */
void EscGadgetWidget::disconnectPort()
{
    if (escSerial)
    {
        delete escSerial;
        escSerial = NULL;
    }

    connectedStatus = false;
}

/**
  * @brief Connect to the selected port
  */
void EscGadgetWidget::connectPort()
{
    int device_idx = m_widget->comboBox_Ports->currentIndex();
    if (device_idx < 0)
        return;

    QString device_str = m_widget->comboBox_Ports->currentText().trimmed();
    Q_ASSERT(!device_str.isEmpty());

    int type = NO_PORT;
    if (device_str.toLower().startsWith("com: "))
    {
        type = SERIAL_PORT;
        device_str.remove(0, 5);
        device_str = device_str.trimmed();

        QString str = getSerialPortDevice(device_str);
        if (str.isEmpty())
            return;

        PortSettings settings;
        settings.BaudRate = BAUD115200;
        settings.DataBits = DATA_8;
        settings.Parity = PAR_NONE;
        settings.StopBits = STOP_1;
        settings.FlowControl = FLOW_OFF;
        settings.Timeout_Millisec = 1000;

        QextSerialPort *serial_dev = new QextSerialPort(str, settings);
        if (!serial_dev)
            return;

        if (!serial_dev->open(QIODevice::ReadWrite))
        {
            delete serial_dev;
            return;
        }

        escSerial = new EscSerial(serial_dev);

        refreshStatus();
        refreshConfiguration();
    } else
        return;


    connectedStatus = true;
}

/**
  * @brief Toggle the connection status
  */
void EscGadgetWidget::connectDisconnect()
{
    if (connectedStatus) {
        disconnectPort();
    } else {
        connectPort();
    }

    m_widget->connectButton->setText(connectedStatus ? "Disconnect" : "Connect");

    if(connectedStatus) {
        connected();
    } else {
        disconnected();
    }
}

/**
  * @brief When UAVObject is updated send to EscSettings via serial
  */
void EscGadgetWidget::sendConfiguration(UAVObject*)
{
    Q_ASSERT(escSerial != NULL);

    EscSettings *escSettings = EscSettings::GetInstance(getObjectManager());
    Q_ASSERT(escSettings != NULL);

    EscSettings::DataFields escSettingsData = escSettings->getData();
    escSerial->setSettings(escSettingsData);
}

/**
  * @brief On connection and peridically get status
  */
void EscGadgetWidget::refreshStatus()
{
    Q_ASSERT(escSerial != NULL);

    EscStatus::DataFields escStatusData = escSerial->getStatus();
    EscStatus *escStatus = EscStatus::GetInstance(getObjectManager());
    escStatus->setData(escStatusData);

    m_widget->dutyCycle->display(escStatusData.DutyCycle);
    m_widget->rpm->display(escStatusData.CurrentSpeed);
    m_widget->battery->display(escStatusData.Battery);
    m_widget->current->display(escStatusData.Current);
    m_widget->setPoint->display(escStatusData.SpeedSetpoint);
}

/**
  * @brief On connection and peridically get configuration
  */
void EscGadgetWidget::refreshConfiguration()
{
    Q_ASSERT(escSerial != NULL);

    EscSettings::DataFields escSettingsData = escSerial->getSettings();
    EscSettings *escSettings = EscSettings::GetInstance(getObjectManager());
    escSettings->setData(escSettingsData);

    m_widget->spinRisingKp->setValue(escSettingsData.RisingKp);
    m_widget->spinFallingKp->setValue(escSettingsData.FallingKp);
    m_widget->spinSortCurrentLimit->setValue(escSettingsData.SoftCurrentLimit);
    m_widget->spinHardCurrentLimit->setValue(escSettingsData.HardCurrentLimit);
    m_widget->spinKv->setValue(escSettingsData.Kv);
    m_widget->closedLoopCheck->setChecked(escSettingsData.Mode == EscSettings::MODE_CLOSED);
    m_widget->spinPhase->setValue(escSettingsData.CommutationPhase);
}

/**
  * @brief Tell the ESC to save the settings to the flash
  */
void EscGadgetWidget::saveConfiguration()
{
    Q_ASSERT(escSerial != NULL);

    applyConfiguration();;
    escSerial->saveSettings();
}

/**
  * @brief Get the settings from the UI and send to the ESC
  */
void EscGadgetWidget::applyConfiguration()
{
    EscSettings *escSettings = EscSettings::GetInstance(getObjectManager());
    Q_ASSERT(escSettings != NULL);

    EscSettings::DataFields escSettingsData = escSettings->getData();
    escSettingsData.RisingKp = m_widget->spinRisingKp->value();
    escSettingsData.FallingKp = m_widget->spinFallingKp->value();
    escSettingsData.SoftCurrentLimit = m_widget->spinSortCurrentLimit->value();
    escSettingsData.HardCurrentLimit = m_widget->spinHardCurrentLimit->value();
    escSettingsData.Kv = m_widget->spinKv->value();
    escSettingsData.Mode = m_widget->closedLoopCheck->checkState() ? EscSettings::MODE_CLOSED : EscSettings::MODE_OPEN;
    escSettings->setData(escSettingsData); // Triggers the update to send
    escSettingsData.CommutationPhase = m_widget->spinPhase->value();
}

/**
  * When connected enable UI
  */
void EscGadgetWidget::connected()
{
    connect(EscSettings::GetInstance(getObjectManager()),
            SIGNAL(objectUpdated(UAVObject*)), this,
            SLOT(sendConfiguration(UAVObject*)));
    connect(&refreshTimer, SIGNAL(timeout()), this, SLOT(refreshStatus()));
    refreshTimer.start(1000);

    m_widget->saveSettingsToFlash->setEnabled(true);
    m_widget->saveSettingsToRAM->setEnabled(true);
}

/**
  * When disconnected disable UI and signals
  */
void EscGadgetWidget::disconnected()
{
    disconnect(EscSettings::GetInstance(getObjectManager()),
            SIGNAL(objectUpdated(UAVObject*)), this,
            SLOT(sendConfiguration(UAVObject*)));
    disconnect(&refreshTimer, SIGNAL(timeout()), this, SLOT(refreshStatus()));
    refreshTimer.stop();

    m_widget->saveSettingsToFlash->setEnabled(false);
    m_widget->saveSettingsToRAM->setEnabled(false);
}

/**
  * @brief Return handle to object manager
  */
UAVObjectManager * EscGadgetWidget::getObjectManager()
{
    ExtensionSystem::PluginManager *pm = ExtensionSystem::PluginManager::instance();
    UAVObjectManager * objMngr = pm->getObject<UAVObjectManager>();
    Q_ASSERT(objMngr);
    return objMngr;
}
