include(../defaults.pri)
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXX = icpc
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

SOURCES += \
    main.cpp \

INCLUDEPATH += $$PWD/../src

LIBS += -L../src/
LIBS += -lmolecular-dynamics

icpc {
    copydata.commands = cp -f $$OUT_PWD/../src/*.so* $$OUT_PWD/
} else {
    copydata.commands = cp -f $$OUT_PWD/../src/*.dylib $$OUT_PWD/
}

first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata
