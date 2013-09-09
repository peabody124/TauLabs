/**
 ******************************************************************************
 * @file       leapcontrolgadget.h
 * @author     Tau Labs, http://taulabs.org, Copyright (C) 2013
 * @addtogroup GCSPlugins GCS Plugins
 * @{
 * @addtogroup LeapControl Leap Controller plugin
 * @{
 * @brief A gadget to control the UAV from a Leap Controller
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

#ifndef LEAPCONTROLGADGET_H
#define LEAPCONTROLGADGET_H

#include <coreplugin/iuavgadget.h>
#include "uavtalk/telemetrymanager.h"
#include "leapcontrolplugin.h"
#include "leapcontrol.h"
#include <QTimer>
#include <QTime>

namespace Core {
class IUAVGadget;
}

class LeapControlGadgetWidget;

using namespace Core;

class LeapControlGadget : public Core::IUAVGadget
{
    Q_OBJECT
public:
    LeapControlGadget(QString classId, LeapControlGadgetWidget *widget, QWidget *parent = 0, QObject *plugin=0);
    ~LeapControlGadget();

    QWidget *widget() { return m_widget; }
    QString contextHelpId() const { return QString(); }

    void loadConfiguration(IUAVGadgetConfiguration* config);

private:
    QWidget     * m_widget;
    LeapControl * m_leapControl;

signals:

private slots:
    void onAutopilotConnect();
    void onAutopilotDisconnect();

protected slots:
    void handUpdated(bool present, double x, double y, double z, double roll, double pitch, double yaw);
};


#endif // LEAPCONTROLGADGET_H

/**
 * @}
 * @}
 */
