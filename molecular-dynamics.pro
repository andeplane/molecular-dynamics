TEMPLATE = subdirs
CONFIG += ordered

CONFIG += c++11

SUBDIRS += \
    src \
    #tests \
    example \
    createMovie
    #Benchmarks

BUILD_DIR = $$OUT_PWD
