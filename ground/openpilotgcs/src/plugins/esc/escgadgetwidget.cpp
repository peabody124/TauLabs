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
    m_ioDevice(NULL),
    connected(false)
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

    connect(m_widget->connectButton, SIGNAL(clicked()), this, SLOT(connectDisconnect()));
//    connect(m_widget->refreshPorts, SIGNAL(clicked()), this, SLOT(getPorts()));
//    connect(m_widget->pushButton_Save, SIGNAL(clicked()), this, SLOT(saveToFlash()));
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
    if (m_ioDevice)
    {
        m_ioDevice->close();
        disconnect(m_ioDevice, 0, 0, 0);
        delete m_ioDevice;
        m_ioDevice = NULL;
    }

    connected = false;
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

        m_ioDevice = serial_dev;

        qDebug() << "Connected";

    } else
        return;


    connected = true;
}

/**
  * @brief Toggle the connection status
  */
void EscGadgetWidget::connectDisconnect()
{
    if (connected) {
        qDebug() << "Disconnecting";
        disconnectPort();
    } else {
        qDebug() << "Connecting";
        connectPort();
    }

    m_widget->connectButton->setText(connected ? "Disconnect" : "Connect");
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
