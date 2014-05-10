include(../defaults.pri)
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += \
    main.cpp \
    forces.cpp \
    system.cpp \
    generator.cpp \
    units.cpp

INCLUDEPATH += /repos/UnitTest++/src
INCLUDEPATH += $$PWD/../src

LIBS += -L/repos/UnitTest++/ -L../src/
LIBS += -lUnitTest++ -lmolecular-dynamics

icpc {
    copydata.commands = cp -f $$OUT_PWD/../src/*.so* $$OUT_PWD/
} else {
    copydata.commands = cp -f $$OUT_PWD/../src/*.dylib $$OUT_PWD/
}

first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata
