TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += \
    main.cpp \
    forces.cpp \
    system.cpp

INCLUDEPATH += /repos/UnitTest++/src
INCLUDEPATH += $$PWD/../src

LIBS += -L/repos/UnitTest++/ -L../src/
LIBS += -lUnitTest++ -lmolecular-dynamics

copydata.commands = cp -f $$OUT_PWD/../src/*.dylib $$OUT_PWD/
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata
