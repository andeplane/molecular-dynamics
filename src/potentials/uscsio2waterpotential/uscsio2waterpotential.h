#pragma once
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <atomiterators/atomiterators.h>
#include <particles/customproperty.h>
#include <potentials/uscsio2waterpotential/uscpotentialparameters.h>
#include <potentials/uscsio2waterpotential/usceffectiveforcefieldinterpolater.h>

typedef vector<double> potential_energy_at_cutoff_t;
typedef vector<vector<vector<int> > > configuration_map_t;

enum class AtomConfiguration {NotUsed = 0, Si_O=1, Si_Si=2, O_O=3, O_H=4, H_H=5, Si_H=6, O_Si_O=7, Si_O_Si=8, H_O_H=9, Si_O_H=10, NumberOfConfigurations=6};
inline int operator + (AtomConfiguration val) {
    return static_cast<int>(val);
}

class USCSIO2WaterPotential : public Potential
{
private:
    AtomIteratorDefault m_iteratorDefault;
    USCWaterPotentialParameters *m_waterParameters;
    USCWaterPotentialParameters *m_silicaParameters;
    // Two particle coefficients
    double m_twoParticleCutoffDistance;
    double m_threeParticleCutoffDistance;

    double m_deltaR2; // Step length in precomputed table
    double m_oneOverDeltaR2; // Step length in precomputed table
    int m_numberOfPrecomputedTwoParticleForces;

    vector<int> m_atomicNumberMap;
    configuration_map_t m_configurationMap;
    vector<int> m_activeAtomTypes;
    void initialize();
    void calculatePrecomputedTwoParticleForces();
    inline void sortAtoms(Atom *&atom1, Atom *&atom2, Atom *&atom3, int atomConfiguration) {
        // We expect atom2 to be the "special" atom in the configuration
        if(atomConfiguration == +AtomConfiguration::O_Si_O) {
            if(atom2->type()->atomicNumber() == +AtomTypes::Silicon) return;
            else if(atom1->type()->atomicNumber() == +AtomTypes::Silicon) std::swap(atom2,atom1);
            else std::swap(atom2,atom3);
        } else {
            if(atom2->type()->atomicNumber() == +AtomTypes::Oxygen) return;
            else if(atom1->type()->atomicNumber() == +AtomTypes::Oxygen) std::swap(atom2,atom1);
            else std::swap(atom2,atom3);
        }
    }
public:
    USCSIO2WaterPotential();
    void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
    void threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3);
    std::string coefficientString() const;
    int numberOfPrecomputedTwoParticleForces() const;
    void setNumberOfPrecomputedTwoParticleForces(int numberOfPrecomputedTwoParticleForces);
};
