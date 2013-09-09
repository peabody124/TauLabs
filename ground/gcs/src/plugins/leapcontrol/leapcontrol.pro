TEMPLATE = lib
TARGET = LeapControl
QT += svg
QT += opengl
QT += network

include(../../taulabsgcsplugin.pri) 
include(../../plugins/coreplugin/coreplugin.pri) 
include(../../plugins/uavobjects/uavobjects.pri)
include(../../plugins/uavtalk/uavtalk.pri)
include(leap.pri)

HEADERS += leapcontrolgadget.h \
    leapinterface.h
HEADERS += leapcontrolgadgetwidget.h
HEADERS += leapcontrolgadgetfactory.h
HEADERS += leapcontrolplugin.h

SOURCES += leapcontrolgadget.cpp \
    leapinterface.cpp
SOURCES += leapcontrolgadgetwidget.cpp
SOURCES += leapcontrolgadgetfactory.cpp
SOURCES += leapcontrolplugin.cpp

OTHER_FILES += LeapControl.pluginspec \

FORMS += leapcontrol.ui
