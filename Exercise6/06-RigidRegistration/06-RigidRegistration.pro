TARGET = RigidRegistration
TEMPLATE = app

CONFIG -= qt
CONFIG -= app_bundle
CONFIG += console


SOURCES += RegistrationViewer.cc \
    main.cc \
    Transformation.cc

HEADERS += RegistrationViewer.hh \
    Transformation.hh \
    Fitting.hh

include(../LibsInclude.pri)
