#pragma once
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <atomiterators/atomiteratordefault.h>
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
typedef vector<double> coeff_r1s;
typedef vector<double> coeff_r4s;
typedef vector<double> coeff_oneOverR1s;
typedef vector<double> coeff_oneOverR4s;
typedef vector<double> coeff_cutoff_distance;

enum class AtomConfiguration {NotUsed = 0, Si_O=1, Si_Si=2, O_O=3, O_Si_O=4, Si_O_Si=5, NumberOfConfigurations=6};
int operator + (AtomConfiguration val);


class USCSIO2Potential : public Potential
{
private:
    friend std::ostream& operator<<(std::ostream&stream, const USCSIO2Potential&atom);
    AtomIteratorDefault m_iteratorDefault;
    // Two particle coefficients
    double m_maxCutoffDistance;
    coeff_cutoff_distance cutoffDistances;
    coeff_cutoff_distance cutoffDistancesSquared;
    coeff_eta_ij eta_ij;
    coeff_H_ij H_ij;
    coeff_D_ij D_ij;
    coeff_W_ij W_ij;
    coeff_Z_i Z_i;

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
    void initialize();
    void sortAtoms(Atom *&atom1, Atom *&atom2, Atom *&atom3, int atomConfiguration);

public:
    USCSIO2Potential();
    virtual void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
    void threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3);
    double maxCutoffDistance() const;
    std::string coefficientString() const;
};
