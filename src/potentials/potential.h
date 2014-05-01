#ifndef POTENTIAL_H
#define POTENTIAL_H

class AtomManager;
class Potential
{
public:
    Potential();
    virtual void calculateForces(AtomManager &atomManager) = 0;
};

#endif // POTENTIAL_H
