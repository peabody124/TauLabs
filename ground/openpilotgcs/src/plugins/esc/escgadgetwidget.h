/**
 ******************************************************************************
 *
 * @file       pipxtremegadgetwidget.h
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

#ifndef ESCGADGETWIDGET_H
#define ESCGADGETWIDGET_H

#include "ui_esc.h"

#include <qextserialport.h>
#include <qextserialenumerator.h>

#include "uavtalk/telemetrymanager.h"
#include "extensionsystem/pluginmanager.h"
#include "uavobjects/uavobjectmanager.h"
#include "uavobjects/uavobject.h"

#include "coreplugin/icore.h"
#include "coreplugin/connectionmanager.h"

#include "rawhid/rawhidplugin.h"

#include <QtGui/QWidget>
#include <QLabel>
#include <QLineEdit>
#include <QThread>
#include <QMessageBox>
#include <QTimer>
#include <QtCore/QVector>
#include <QtCore/QIODevice>
#include <QtCore/QLinkedList>
#include <QMutex>
#include <QMutexLocker>
#include "widgetbar.h"
#include "escserial.h"

class EscGadgetWidget : public QWidget
{
    Q_OBJECT

public:
    EscGadgetWidget(QWidget *parent = 0);
    ~EscGadgetWidget();

public slots:
    void onComboBoxPorts_currentIndexChanged(int index);

protected:
    void resizeEvent(QResizeEvent *event);

private:
    Ui_EscWidget	*m_widget;

    EscSerial       *escSerial;

    QString getSerialPortDevice(const QString &friendName);

    QTimer refreshTimer;

    void getSettings();
    void sendSettings();

    UAVObjectManager * getObjectManager();

    bool connectedStatus;

private slots:
    void connectDisconnect();
    void disconnectPort();
    void connectPort();
    void getPorts();

    void saveConfiguration();
    void applyConfiguration();
    void sendConfiguration(UAVObject*);
    void refreshStatus();
    void refreshConfiguration();

    void connected();
    void disconnected();

    void queryDfuDevice();
    void rescueCode();
    void updateCode();
};

#endif
