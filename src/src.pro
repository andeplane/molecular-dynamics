include(../defaults.pri)
TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
TARGET = molecular-dynamics

HEADERS += \
    particles/atom.h \
    filemanager/filemanager.h \
    generator.h \
    simulator.h \
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
    utils/utils.h \
    statistics/statisticalproperty.h \
    statistics/kineticenergysampler.h \
    statistics/statisticalvalue.h \
    statistics/temperaturesampler.h \
    statistics/potentialenergysampler.h \
    node.h \
    statistics/totalenergysampler.h \
    statistics/paircorrelationsampler.h \
    statistics/neighborlist.h \
    statistics/statistics.h \
    outputs/standardconsoleoutput.h \
    particles/particle.h \
    modifiers/modifier.h \
    modifiers/modifiers.h \
    modifiers/berendsenthermostat.h \
    particles/customproperty.h \
    utils/vec3.h \
    potentials/uscsio2waterpotential/uscsio2waterpotential.h \
    potentials/uscsio2waterpotential/uscpotentialparameters.h \
    potentials/uscsio2waterpotential/usceffectiveforcefieldinterpolater.h

SOURCES += \
    particles/atom.cpp \
    filemanager/filemanager.cpp \
    generator.cpp \
    simulator.cpp \
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
    utils/utils.cpp \
    statistics/statisticalproperty.cpp \
    statistics/kineticenergysampler.cpp \
    statistics/temperaturesampler.cpp \
    statistics/potentialenergysampler.cpp \
    node.cpp \
    statistics/totalenergysampler.cpp \
    statistics/paircorrelationsampler.cpp \
    statistics/neighborlist.cpp \
    outputs/standardconsoleoutput.cpp \
    particles/particle.cpp \
    modifiers/modifier.cpp \
    modifiers/berendsenthermostat.cpp \
    particles/customproperty.cpp \
    utils/vec3.cpp \
    potentials/uscsio2waterpotential/uscsio2waterpotential.cpp \
    potentials/uscsio2waterpotential/uscpotentialparameters.cpp \
    potentials/uscsio2waterpotential/usceffectiveforcefieldinterpolater.cpp

icpc {
    QMAKE_LFLAGS += -staticlib
    QMAKE_LFLAGS_SONAME -= -Wl,-soname,
}

# message($$QMAKE_LFLAGS_SONAME)
