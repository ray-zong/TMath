#-------------------------------------------------
#
# Project created by QtCreator 2016-05-31T12:29:03
#
#-------------------------------------------------

QT       -= gui

TARGET = TMath
TEMPLATE = lib

DEFINES += TMATH_LIBRARY

INCLUDEPATH += include

SOURCES += src/TMath.cpp \
    src/TVector.cpp \
    src/TMatrix.cpp

HEADERS += include/TMath.h \
    include/TVector.h \
    include/TMatrix.h \
    include/TGlobal.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
