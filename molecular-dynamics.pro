TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    system.cpp \
    topology.cpp \
    atom.cpp \
    filemanager.cpp \
    potential.cpp \
    lennardjonespotential.cpp \
    integrator.cpp \
    velocityverlet.cpp \
    simulator.cpp \
    statisticssampler.cpp \
    statisticalproperty.cpp \
    unitconverter.cpp \
    generator.cpp \
    systemcell.cpp

HEADERS += \
    system.h \
    topology.h \
    atom.h \
    filemanager.h \
    potential.h \
    lennardjonespotential.h \
    integrator.h \
    velocityverlet.h \
    simulator.h \
    statisticssampler.h \
    statisticalproperty.h \
    unitconverter.h \
    generator.h \
    systemcell.h

