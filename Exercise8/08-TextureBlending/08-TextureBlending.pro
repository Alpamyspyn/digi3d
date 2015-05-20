TARGET = TextureBlending
TEMPLATE = app

CONFIG += qt
QT = core xml gui
CONFIG -= app_bundle
CONFIG += console


SOURCES += TextureBlending.cc \
    VertexVisibility.cc

HEADERS += VertexVisibility.hh

include(../LibsInclude.pri)
