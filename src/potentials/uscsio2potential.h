#pragma once
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <atomiterators/atomiterators.h>
#include <particles/customproperty.h>

typedef vector<double> coeff_B_ijk;
typedef vector<double> coeff_theta_zero;
typedef vector<double> coeff_r0;
typedef vector<double> coeff_ksi;
typedef vector<double> coeff_C_ijk;
typedef vector<int> coeff_eta_ij;
typedef vector<double> coeff_H_ij;
typedef vector<double> coeff_D_ij;
typedef vector<double> coeff_W_ij;
typedef vector<double> coeff_Z_i;
typedef vector<double> coeff_alpha_i;
typedef vector<double> coeff_r1s;
typedef vector<double> coeff_r4s;
typedef vector<double> coeff_oneOverR1s;
typedef vector<double> coeff_oneOverR4s;
typedef vector<double> coeff_cutoff_distance;
typedef vector<double> potential_energy_at_cutoff;
typedef vector<vector<vector<double> > > precomputed_two_particle_forces;
typedef vector<vector<vector<double> > > precomputed_two_particle_potential;


class USCSIO2Potential : public Potential
{
private:
    friend std::ostream& operator<<(std::ostream&stream, const USCSIO2Potential &potential);
    AtomIteratorDefault m_iteratorDefault;
    AtomIteratorAllPairs m_iteratorAllPairs;

    // Two particle coefficients
    double m_maxTwoParticleCutoffDistance;
    double m_maxThreeParticleCutoffDistance;
    coeff_cutoff_distance cutoffDistances;
    coeff_cutoff_distance cutoffDistancesSquared;
    coeff_cutoff_distance oneOverCutoffDistancesSquared;
    coeff_eta_ij eta_ij;
    coeff_H_ij H_ij;
    coeff_D_ij D_ij;
    coeff_W_ij W_ij;
    coeff_Z_i Z_i;
    coeff_alpha_i alpha_i;

    int m_numberOfPrecomputedTwoParticleForces;
    precomputed_two_particle_forces m_precomputedTwoParticleForces;
    precomputed_two_particle_potential m_precomputedTwoParticlePotential;

    potential_energy_at_cutoff m_potentialEnergyAtCutoff; // To correct for cutoff distance
    double m_deltaR2; // Step length in precomputed table
    double m_oneOverDeltaR2; // Step length in precomputed table

    coeff_r1s r1s;
    coeff_r4s r4s;

    coeff_oneOverR1s oneOverR1s;
    coeff_oneOverR4s oneOverR4s;

    // Three particle coefficients
    coeff_B_ijk B_ijk;
    coeff_theta_zero cosThetaZero;
    coeff_r0 r0;
    coeff_ksi ksi;
    coeff_C_ijk C_ijk;

    vector<int> m_atomicNumberMap;
    vector<vector<vector<int> > > m_configurationMap;
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
    USCSIO2Potential();
    virtual void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
    void threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3);
    std::string coefficientString() const;
    int numberOfPrecomputedTwoParticleForces() const;
    void setNumberOfPrecomputedTwoParticleForces(int numberOfPrecomputedTwoParticleForces);
};
