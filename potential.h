#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <vector>

class Atom;
using std::vector;

class Potential
{
public:
    Potential();
    virtual void calculateForces(vector<Atom> &atoms) = 0;
};

#endif // POTENTIAL_H
