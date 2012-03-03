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
	widgetbar.h \
    escserial.h \
    stm32.h

SOURCES += escgadget.cpp \
    escgadgetconfiguration.cpp \
    escgadgetfactory.cpp \
    escgadgetoptionspage.cpp \
    escgadgetwidget.cpp \
    escplugin.cpp \
	widgetbar.cpp \
    escserial.cpp \
    stm32.cpp

OTHER_FILES += ESC.pluginspec

FORMS += \
    esc.ui \
    escgadgetoptionspage.ui

RESOURCES += \
    esc.qrc




