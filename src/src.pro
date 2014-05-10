include(../defaults.pri)
TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
TARGET = molecular-dynamics

HEADERS += \
    atom.h \
    filemanager/filemanager.h \
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
    includes.h \
    random.h \
    integrators/integrators.h \
    atomiterators/atomiterators.h \
    potentials/potentials.h \
    atomiterators/atomiteratordefault.h \
    atomiterators/atomiteratorallpairs.h \
    potentials/uscsio2potential.h \
    filemanager/mts0io.h \
    utils/utils.h

SOURCES += \
    atom.cpp \
    filemanager/filemanager.cpp \
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
    atomlist.cpp \
    random.cpp \
    atomiterators/atomiteratordefault.cpp \
    atomiterators/atomiteratorallpairs.cpp \
    potentials/uscsio2potential.cpp \
    filemanager/mts0io.cpp \
    utils/utils.cpp

icpc {
    QMAKE_LFLAGS += -staticlib
    QMAKE_LFLAGS_SONAME -= -Wl,-soname,
}

message($$QMAKE_LFLAGS_SONAME)
