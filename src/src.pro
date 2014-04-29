TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

TARGET = molecular-dynamics

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
    potentials/potential.h \
    atomtype.h \
    atomlist.h

SOURCES += \
    atom.cpp \
    filemanager.cpp \
    generator.cpp \
    main.cpp \
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
    potentials/potential.cpp \
    atomtype.cpp \
    atomlist.cpp
