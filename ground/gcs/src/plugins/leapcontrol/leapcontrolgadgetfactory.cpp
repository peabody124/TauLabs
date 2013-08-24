/**
 ******************************************************************************
 * @file       leapcontrolgadgetfactory.cpp
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

#include "leapcontrolgadgetfactory.h"
#include "leapcontrolgadgetwidget.h"
#include "leapcontrolgadget.h"
#include <coreplugin/iuavgadget.h>

LeapControlGadgetFactory::LeapControlGadgetFactory(QObject *parent) :
        IUAVGadgetFactory(QString("LeapControlGadget"),
                          tr("Leap"),
                          parent)
{
}

LeapControlGadgetFactory::~LeapControlGadgetFactory()
{

}

IUAVGadget* LeapControlGadgetFactory::createGadget(QWidget *parent) {
    LeapControlGadgetWidget* gadgetWidget = new LeapControlGadgetWidget(parent);
    return new LeapControlGadget(QString("LeapControlGadget"), gadgetWidget, parent, this->parent());
}

/**
 * @}
 * @}
 */

