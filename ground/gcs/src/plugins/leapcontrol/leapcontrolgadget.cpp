/**
 ******************************************************************************
 * @file       leapcontrolgadget.cpp
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
#include "leapcontrolgadget.h"
#include "leapcontrolgadgetwidget.h"
#include "leapinterface.h"

#include "extensionsystem/pluginmanager.h"
#include "uavobjectmanager.h"
#include "uavobject.h"
#include <QDebug>

LeapControlGadget::LeapControlGadget(QString classId, LeapControlGadgetWidget *widget, QWidget *parent, QObject *plugin) :
        IUAVGadget(classId, parent),
        m_widget(widget)
{
    Q_UNUSED(plugin)

    ExtensionSystem::PluginManager *pm = ExtensionSystem::PluginManager::instance();
    Q_ASSERT(pm != NULL);
    LeapInterface *lf = pm->getObject<LeapInterface>();
    Q_ASSERT(lf != NULL);

    UAVObjectManager *objmngr = pm->getObject<UAVObjectManager>();
    Q_ASSERT(objmngr != NULL);

    m_leapControl = LeapControl::GetInstance(objmngr);
    Q_ASSERT(m_leapControl);

    // Listen to autopilot connection events
    TelemetryManager* telMngr = pm->getObject<TelemetryManager>();
    connect(telMngr, SIGNAL(connected()), this, SLOT(onAutopilotConnect()));
    connect(telMngr, SIGNAL(disconnected()), this, SLOT(onAutopilotDisconnect()));

    connect(lf, SIGNAL(handUpdated(bool,double,double,double,double,double,double)), this,
            SLOT(handUpdated(bool,double,double,double,double,double,double)));
}

LeapControlGadget::~LeapControlGadget()
{
    delete m_widget;
}

void LeapControlGadget::loadConfiguration(IUAVGadgetConfiguration* config)
{
    Q_UNUSED(config);
}

void LeapControlGadget::onAutopilotConnect()
{
    // By default, the LeapControl object metadata does not set it up to send updates
    // over the telemetry link, so we gotta update this metadata when we instanciate the gadget.
    // A 30ms update rate requires a hefty telemetry link, but that's what you should have if
    // you are playing with this :
    UAVObject::Metadata mdata;
    mdata = m_leapControl->getMetadata();

    UAVObject::SetGcsTelemetryUpdateMode(mdata, UAVObject::UPDATEMODE_PERIODIC);
    mdata.gcsTelemetryUpdatePeriod = 30;

    UAVObject::SetFlightTelemetryUpdateMode(mdata, UAVObject::UPDATEMODE_PERIODIC);
    mdata.flightTelemetryUpdatePeriod = 0;

    m_leapControl->setMetadata(mdata);

}

void LeapControlGadget::onAutopilotDisconnect()
{
    // We want to reset metadata to default values once we
    // disconnect
    UAVObject::Metadata mdata;
    mdata = m_leapControl->getDefaultMetadata();
    m_leapControl->setMetadata(mdata);
}

//! When hand is updated, store the information in a UAVO
void LeapControlGadget::handUpdated(bool present, double /*x*/, double /*y*/, double /*z*/,
                                    double roll, double pitch, double yaw)
{
    LeapControl::DataFields leapControlData = m_leapControl->getData();

    if (present) {
        leapControlData.Present = LeapControl::PRESENT_TRUE;
        leapControlData.Roll = roll;
        leapControlData.Pitch = pitch;
        leapControlData.Yaw = yaw;
    } else {
        leapControlData.Present = LeapControl::PRESENT_FALSE;
    }

    m_leapControl->setData(leapControlData);
}

/**
 * @}
 * @}
 */
