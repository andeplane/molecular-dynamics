TEMPLATE = subdirs
CONFIG += ordered

CONFIG += icpc

SUBDIRS += \
    src \
    #tests \
    example \
    createMovie \
    #Benchmarks

BUILD_DIR = $$OUT_PWD
