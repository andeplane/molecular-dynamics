#ifndef LENNARDJONESPOTENTIAL_H
#define LENNARDJONESPOTENTIAL_H
#include <potentials/potential.h>
#include <atomiterators/atomiteratordefault.h>
#include <atomiterators/atomiteratorallpairs.h>

class LennardJonesPotential : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
    double m_cutoffDistance;
    bool m_calculateForcesBetweenAllPairsWithMinimumImageConvention;
    unsigned long m_numberOfComputedPairs;
    unsigned long m_numberOfComputedPairsWithinCutoffDistance;
    double m_potentialEnergyCorrection; // To correct for cutoff distance
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
    void twoParticleActionMinimumImageConvention(Atom *atom1, Atom *atom2);
    unsigned long numberOfComputedPairs() const;
    unsigned long numberOfComputedPairsWithinCutoffDistance() const;
    bool calculateForcesBetweenAllPairsWithMinimumImageConvention() const;
    void setCalculateForcesBetweenAllPairsWithMinimumImageConvention(bool calculateForcesBetweenAllPairsWithMinimumImageConvention);
};

#endif // LENNARDJONESPOTENTIAL_H
