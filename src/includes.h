#pragma once

#include <system.h>
#include <simulator.h>
#include <potentials/potential.h>
#include <integrators/integrator.h>
#include <atom.h>
#include <generator.h>
#include <atomlist.h>
#include <atommanager.h>
#include <potentials/potentials.h>
#include <integrators/integrators.h>
#include <atomiterators/atomiterators.h>
#include <unitconverter.h>
#include <vector>

template<typename T>
std::ostream& operator<<(std::ostream &stream, const std::vector<T> &vec) {
    stream << "[";
    for(unsigned long i=0; i<vec.size(); i++) {
        if(i+1 == vec.size()) stream << vec[i];
        else stream << vec[i] << ", ";
    }
    stream << "]";

    return stream;
}

// #define DEBUG
template <typename T>
T &safeOrQuickVectorLookup(vector<T> &vec, int index) {
#ifdef DEBUG
    return vec.at(index);
#else
    return vec[index];
#endif
}
