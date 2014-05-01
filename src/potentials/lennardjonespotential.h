#ifndef LENNARDJONESPOTENTIAL_H
#define LENNARDJONESPOTENTIAL_H
#include <potentials/potential.h>
#include <atomiterators/atomiteratordefault.h>

class LennardJonesPotential : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
    double m_cutoffDistance;
    AtomIteratorDefault m_iterator;
public:
    LennardJonesPotential();
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
    double cutoffDistance() const;
    void setCutoffDistance(double cutoffDistance);
    virtual void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
};

#endif // LENNARDJONESPOTENTIAL_H
