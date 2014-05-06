#pragma once
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <atomiterators/atomiteratordefault.h>
typedef double coeff_eta_ij;

class USCSIO2Potential : public Potential
{
private:
    vector<double> m_cutoffDistances;
    vector<double> m_cutoffDistancesSquared;
    AtomIteratorDefault m_iteratorDefault;
    double m_maxCutoffDistance;
    vector<vector<coeff_eta_ij> > eta_ij;
    vector<int> m_atomicNumberMap;

    void initialize();
public:
    USCSIO2Potential();
    virtual void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
    void threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3);
    vector<double> cutoffDistances() const;
    double maxCutoffDistance() const;
};
