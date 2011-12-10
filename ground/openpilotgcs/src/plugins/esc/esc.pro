TEMPLATE = lib
TARGET = ESC

QT += svg
QT += opengl

include(../../openpilotgcsplugin.pri)
include(../../plugins/coreplugin/coreplugin.pri)
include(../../plugins/uavobjects/uavobjects.pri)
include(../../plugins/uavtalk/uavtalk.pri)
include(../../plugins/rawhid/rawhid.pri)

INCLUDEPATH += ../../libs/qextserialport/src

HEADERS += escgadget.h \
    escgadgetconfiguration.h \
    escgadgetfactory.h \
    escgadgetoptionspage.h \
    escgadgetwidget.h \
    escplugin.h \
	widgetbar.h

SOURCES += escgadget.cpp \
    escgadgetconfiguration.cpp \
    escgadgetfactory.cpp \
    escgadgetoptionspage.cpp \
    escgadgetwidget.cpp \
    escplugin.cpp \
	widgetbar.cpp

OTHER_FILES += esc.pluginspec

FORMS += \
    esc.ui \
    escgadgetoptionspage.ui

RESOURCES += \
    esc.qrc
