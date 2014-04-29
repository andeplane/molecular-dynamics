#ifndef POTENTIAL_H
#define POTENTIAL_H
class AtomManager;

enum class PotentialType { LennardJones };

class Potential
{
public:
    Potential();
    virtual void calculateForces(AtomManager &atomList) = 0;
};

#endif // POTENTIAL_H
