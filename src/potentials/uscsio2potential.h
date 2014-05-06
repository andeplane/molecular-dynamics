#pragma once
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <atomiterators/atomiteratordefault.h>
typedef vector<vector<vector<double> > > coeff_B_ijk;
typedef vector<vector<vector<double> > > coeff_theta_zero;
typedef vector<vector<vector<double> > > coeff_r_0;
typedef vector<vector<vector<double> > > coeff_ksi;
typedef vector<vector<vector<double> > > coeff_C_ijk;
typedef vector<vector<double> > coeff_eta_ij;
typedef vector<vector<double> > coeff_H_ij;
typedef vector<vector<double> > coeff_D_ij;
typedef vector<vector<double> > coeff_W_ij;
typedef vector<vector<double> > coeff_cutoff_distance;
typedef vector<double>          coeff_Z_i;
typedef vector<double>          coeff_r_1s;
typedef vector<double>          coeff_r_4s;

class USCSIO2Potential : public Potential
{
private:
    AtomIteratorDefault m_iteratorDefault;
    // Two particle coefficients
    double m_maxCutoffDistance;
    coeff_cutoff_distance cutoffDistances;
    coeff_eta_ij eta_ij;
    coeff_H_ij H_ij;
    coeff_D_ij D_ij;
    coeff_W_ij W_ij;
    coeff_Z_i Z_i;
    coeff_r_1s r1s;
    coeff_r_4s r4s;

    // Three particle coefficients
    coeff_B_ijk B_ijk;
    coeff_theta_zero theta_zero;
    coeff_r_0 r_0;
    coeff_ksi ksi;
    coeff_C_ijk C_ijk;

    vector<int> m_atomicNumberMap;

    void initialize();
public:
    USCSIO2Potential();
    virtual void calculateForces(AtomManager &atomManager);
    void twoParticleAction(Atom *atom1, Atom *atom2);
    void threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3);
    double maxCutoffDistance() const;
};
