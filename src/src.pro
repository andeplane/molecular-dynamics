TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    atom.cpp \
    filemanager.cpp \
    generator.cpp \
    simulator.cpp \
    statisticalproperty.cpp \
    statisticssampler.cpp \
    system.cpp \
    systemcell.cpp \
    topology.cpp \
    unitconverter.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    potentials/lennardjonespotential.cpp \
    potentials/potential.cpp

HEADERS += \
    atom.h \
    filemanager.h \
    generator.h \
    simulator.h \
    statisticalproperty.h \
    statisticssampler.h \
    system.h \
    systemcell.h \
    topology.h \
    unitconverter.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    potentials/lennardjonespotential.h \
    potentials/potential.h

