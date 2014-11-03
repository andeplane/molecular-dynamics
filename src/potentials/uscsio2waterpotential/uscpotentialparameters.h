#pragma once

#include <vector>
using std::vector;

typedef vector<double> coeff_B_ijk;
typedef vector<double> coeff_theta_zero;
typedef vector<double> coeff_r0;
typedef vector<double> coeff_ksi;
typedef vector<int> coeff_eta_ij;
typedef vector<double> coeff_H_ij;
typedef vector<double> coeff_D_ij;
typedef vector<double> coeff_W_ij;
typedef vector<double> coeff_Z_i;
typedef vector<double> coeff_r4s;
typedef vector<vector<vector<double> > > precomputed_two_particle_forces_t;
typedef vector<vector<vector<double> > > precomputed_two_particle_potential_t;
typedef vector<vector<vector<int> > > configuration_map_t;

enum class PotentialClass {Silica=0, Water=1};
enum class ForceOrPotential {Potential=0, Force=1};

class USCWaterPotentialParameters
{
private:
    precomputed_two_particle_forces_t m_precomputedTwoParticleForces;
    precomputed_two_particle_potential_t m_precomputedTwoParticlePotential;

    vector<int> m_activeAtomTypes;
    vector<int> m_atomicNumberMap;
    configuration_map_t m_configurationMap;

    void setSilicaParameters();
    void setCommonParameters();
    void setWaterParameters();
    double stericRepulsion(ForceOrPotential forceOrPotential, int atomConfiguration, double r);
    double coulomb(ForceOrPotential forceOrPotential, int atomConfiguration, int atomicNumber1, int atomicNumber2, double r);
    double chargeDipole(ForceOrPotential forceOrPotential, int atomConfiguration, double r);
public:
    USCWaterPotentialParameters(PotentialClass c, vector<int> &activeAtomTypes, vector<int> &atomicNumberMap, vector<vector<vector<int> > > &configurationMap);

    // Two particle coefficients
    double cutoffDistance;
    coeff_eta_ij eta_ij;
    coeff_H_ij H_ij;
    coeff_D_ij D_ij;
    coeff_W_ij W_ij;
    coeff_Z_i Z_i;

    double r1s;
    coeff_r4s r4s;

    // Three particle coefficients
    coeff_B_ijk B_ijk;
    coeff_theta_zero cosThetaZero;
    coeff_r0 r0;
    coeff_ksi ksi;

    vector<double> &forcesTable(int atomicNumber1, int atomicNumber2);
    vector<double> &potentialTable(int atomicNumber1, int atomicNumber2);
    void createPrecomputedTables(int numberOfPrecomputedForces, double &deltaR2);
};
