TARGET = RigidRegistration
TEMPLATE = app

CONFIG += qt
QT = core network
CONFIG -= app_bundle
CONFIG += console #c++11


SOURCES += RegistrationViewer.cc \
    Transformation.cc \
    fs_tcp_client.cpp \
    main.cpp \
    fsbinarystream.cpp

HEADERS += RegistrationViewer.hh \
    Transformation.hh \
    Fitting.hh \
    animation_loader.h \
    fs_tcp_client.h \
    fsbinarystream.h

include(../LibsInclude.pri)
