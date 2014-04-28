#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <vector>
using std::vector;

class System;

enum class PotentialType { LennardJones };

class Potential
{
public:
    Potential();
    virtual void calculateForces(System &system) = 0;
};

#endif // POTENTIAL_H
