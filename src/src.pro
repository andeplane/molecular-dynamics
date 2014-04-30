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
    topology.h \
    unitconverter.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    potentials/lennardjonespotential.h \
    potentials/potential.h \
    atomtype.h \
    atomiterators/atomiterator.h \
    atommanager.h \
    cell.h \
    atomlist.h \
    includes.h

SOURCES += \
    atom.cpp \
    filemanager.cpp \
    generator.cpp \
    simulator.cpp \
    statisticalproperty.cpp \
    statisticssampler.cpp \
    system.cpp \
    topology.cpp \
    unitconverter.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    potentials/lennardjonespotential.cpp \
    potentials/potential.cpp \
    atomtype.cpp \
    atomiterators/atomiterator.cpp \
    atommanager.cpp \
    cell.cpp \
    atomlist.cpp
