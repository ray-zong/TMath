#-------------------------------------------------
#
# Project created by QtCreator 2016-05-31T12:29:03
#
#-------------------------------------------------

QT       -= gui

CONFIG(debug, debug|release){
     mac: TARGET = $$join(TARGET,,,_debug)
     win32: TARGET = $$join(TARGET,,,d)
}

TEMPLATE = lib

DEFINES += TMATH_LIBRARY

INCLUDEPATH += include

SOURCES += src/TMath.cpp

HEADERS += include/TMath.h \
    include/TVector.h \
    include/TMatrix.h \
    include/TGlobal.h \
    include/TConstants.h \
    include/TQuaternion.h \
    include/TUtilities.h \
    include/TGLSL_mappings.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
