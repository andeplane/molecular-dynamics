TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    topology.cpp \
    atom.cpp \
    filemanager.cpp \
    potential.cpp

HEADERS += \
    system.h \
    topology.h \
    atom.h \
    filemanager.h \
    potential.h

