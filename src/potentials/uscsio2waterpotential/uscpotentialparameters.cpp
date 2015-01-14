#include <potentials/uscsio2waterpotential/uscpotentialparameters.h>
#include <unitconverter.h>
#include <utils/utils.h>
#include <atomtype.h>
#include <cmath>

using CompPhys::Utils::at;

typedef UnitConverter UC;
enum class AtomConfiguration {NotUsed = 0, Si_O=1, Si_Si=2, O_O=3, O_H=4, H_H=5, Si_H=6, O_Si_O=7, Si_O_Si=8, H_O_H=9, Si_O_H=10, NumberOfConfigurations=11};
inline int operator + (AtomConfiguration val) {
    return static_cast<int>(val);
}

USCWaterPotentialParameters::USCWaterPotentialParameters(PotentialClass c, vector<int> &activeAtomTypes, vector<int> &atomicNumberMap, vector<vector<vector<int> > > &configurationMap) :
    m_activeAtomTypes(activeAtomTypes),
    m_atomicNumberMap(atomicNumberMap),
    m_configurationMap(configurationMap)
{
    B_ijk.resize(+AtomConfiguration::NumberOfConfigurations,0);
    cosThetaZero.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r0.resize(+AtomConfiguration::NumberOfConfigurations,0);
    ksi.resize(+AtomConfiguration::NumberOfConfigurations,0);
    eta_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    H_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    D_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    W_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r4s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    Z_i.resize(atomicNumberMap.size(),0);

    if(c == PotentialClass::Silica) setSilicaParameters();
    else setWaterParameters();
    setCommonParameters();
}

void USCWaterPotentialParameters::setCommonParameters() {
    r1s = UC::lengthFromAngstroms(4.43);
    r4s.at(+AtomConfiguration::Si_O) = UC::lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::Si_Si) = UC::lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::O_O) = UC::lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::H_H) = UC::lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::O_H) = UC::lengthFromAngstroms(1.51113);

    cutoffDistance = UC::lengthFromAngstroms(5.5);

    Z_i.at(+AtomTypes::Silicon) = 1.2; // Given in atomic units
    Z_i.at(+AtomTypes::Hydrogen) = 0.32983; // Given in atomic units

    eta_ij.at(+AtomConfiguration::Si_O) = 9;
    eta_ij.at(+AtomConfiguration::Si_Si) = 11;
    eta_ij.at(+AtomConfiguration::H_H) = 9;
    eta_ij.at(+AtomConfiguration::O_H) = 9;

    D_ij.at(+AtomConfiguration::Si_Si) = 0.0;
    D_ij.at(+AtomConfiguration::Si_O) = 3.456*pow(UC::lengthFromAngstroms(1.0),3); // Wrong units in the paper
    D_ij.at(+AtomConfiguration::O_H) = 0.2611*pow(UC::lengthFromAngstroms(1.0),3); // Wrong units in the paper

    H_ij.at(+AtomConfiguration::Si_O) = 78.3143*UC::energyFromEv(1.0)*pow(UC::lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_O));
    H_ij.at(+AtomConfiguration::Si_Si) = 0.39246*UC::energyFromEv(1.0)*pow(UC::lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_Si));
    H_ij.at(+AtomConfiguration::O_H) = 0.61437*UC::energyFromEv(1.0)*pow(UC::lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_H));

    B_ijk.at(+AtomConfiguration::O_Si_O) = UC::energyFromEv(4.993);
    B_ijk.at(+AtomConfiguration::Si_O_Si) = UC::energyFromEv(19.972);
    B_ijk.at(+AtomConfiguration::H_O_H) = UC::energyFromEv(52.9333);
    B_ijk.at(+AtomConfiguration::Si_O_H) = UC::energyFromEv(36.45265);

    cosThetaZero.at(+AtomConfiguration::O_Si_O) = cos(UC::degreesToRadians(109.47));
    cosThetaZero.at(+AtomConfiguration::Si_O_Si) = cos(UC::degreesToRadians(141.0));
    cosThetaZero.at(+AtomConfiguration::H_O_H) = cos(UC::degreesToRadians(97.9476));
    cosThetaZero.at(+AtomConfiguration::Si_O_H) = cos(UC::degreesToRadians(110));

    r0.at(+AtomConfiguration::O_Si_O) = UC::lengthFromAngstroms(2.60);
    r0.at(+AtomConfiguration::Si_O_Si) = UC::lengthFromAngstroms(2.60);
    r0.at(+AtomConfiguration::Si_O_H) = UC::lengthFromAngstroms(2.6);
    r0.at(+AtomConfiguration::H_O_H) = UC::lengthFromAngstroms(1.4);

    ksi.at(+AtomConfiguration::O_Si_O) = UC::lengthFromAngstroms(1.0);
    ksi.at(+AtomConfiguration::Si_O_Si) = UC::lengthFromAngstroms(1.0);
    ksi.at(+AtomConfiguration::H_O_H) = UC::lengthFromAngstroms(0.75);
    ksi.at(+AtomConfiguration::Si_O_H) = UC::lengthFromAngstroms(0.75);
}

vector<double> &USCWaterPotentialParameters::forcesTable(int atomicNumber1, int atomicNumber2)
{
    int mappedAtomicNumber1 = at(m_atomicNumberMap, atomicNumber1);
    int mappedAtomicNumber2 = at(m_atomicNumberMap, atomicNumber2);

    return at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1), mappedAtomicNumber2);
}

vector<double> &USCWaterPotentialParameters::potentialTable(int atomicNumber1, int atomicNumber2)
{
    int mappedAtomicNumber1 = at(m_atomicNumberMap, atomicNumber1);
    int mappedAtomicNumber2 = at(m_atomicNumberMap, atomicNumber2);

    return at(at(m_precomputedTwoParticlePotential,mappedAtomicNumber1), mappedAtomicNumber2);
}

void USCWaterPotentialParameters::setSilicaParameters() {
    eta_ij.at(+AtomConfiguration::O_O) = 7;
    Z_i.at(+AtomTypes::Oxygen) = -0.6;
    D_ij.at(+AtomConfiguration::O_O) = 1.728*pow(UC::lengthFromAngstroms(1.0),3); // Wrong units in the paper
    H_ij.at(+AtomConfiguration::O_O) = 355.5263*UC::energyFromEv(1.0)*pow(UC::lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O));
}

void USCWaterPotentialParameters::setWaterParameters() {
    eta_ij.at(+AtomConfiguration::O_O) = 9;
    Z_i.at(+AtomTypes::Oxygen) = -0.65966;
    H_ij.at(+AtomConfiguration::O_O) = 1965.88*UC::energyFromEv(1.0)*pow(UC::lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O));
    D_ij.at(+AtomConfiguration::O_O) = 2.0887*pow(UC::lengthFromAngstroms(1.0),3); // Wrong units in the paper
}

#include <iostream>
using namespace std;

void USCWaterPotentialParameters::createPrecomputedTables(int numberOfPrecomputedForces, double &deltaR2) {
    m_precomputedTwoParticleForces.resize(m_activeAtomTypes.size(),
                                          vector<vector<double > >(m_activeAtomTypes.size(),
                                                                   vector<double >(numberOfPrecomputedForces+1,0))); // +1 because we need rom for the extra value to interpolate in the last bin

    m_precomputedTwoParticlePotential.resize(m_activeAtomTypes.size(),
                                             vector<vector<double > >(m_activeAtomTypes.size(),
                                                                      vector<double >(numberOfPrecomputedForces+1,0))); // +1 because we need rom for the extra value to interpolate in the last bin


    double rMinSquared = 0;
    double rMaxSquared = cutoffDistance*cutoffDistance;
    deltaR2 = (rMaxSquared - rMinSquared) / (numberOfPrecomputedForces-1);

    for(int atomicNumber1 : m_activeAtomTypes) {
        int mappedAtomicNumber1 = m_atomicNumberMap.at(atomicNumber1);
        for(int atomicNumber2 : m_activeAtomTypes) {
            int mappedAtomicNumber2 = m_atomicNumberMap.at(atomicNumber2);
            int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),+AtomTypes::NoAtom);
            if(atomConfiguration == +AtomConfiguration::NotUsed) continue; // This configuration has zero contribution to the force

            for(int i=0; i<numberOfPrecomputedForces; i++) {
                double r2 = rMinSquared + i*deltaR2;

                if(r2 > cutoffDistance) continue;
                double r = sqrt(r2);

                double force = stericRepulsion(ForceOrPotential::Force, atomConfiguration, r)
                        + coulomb(ForceOrPotential::Force, atomConfiguration, atomicNumber1, atomicNumber2, r)
                        + chargeDipole(ForceOrPotential::Force, atomConfiguration, r);

                double forceAtRCut = stericRepulsion(ForceOrPotential::Force, atomConfiguration, cutoffDistance)
                        + coulomb(ForceOrPotential::Force, atomConfiguration, atomicNumber1, atomicNumber2, cutoffDistance)
                        + chargeDipole(ForceOrPotential::Force, atomConfiguration, cutoffDistance);

                double potential = stericRepulsion(ForceOrPotential::Potential, atomConfiguration, r)
                        + coulomb(ForceOrPotential::Potential, atomConfiguration, atomicNumber1, atomicNumber2, r)
                        + chargeDipole(ForceOrPotential::Potential, atomConfiguration, r);

                double potentialAtRCut = stericRepulsion(ForceOrPotential::Potential, atomConfiguration, cutoffDistance)
                        + coulomb(ForceOrPotential::Potential, atomConfiguration, atomicNumber1, atomicNumber2, cutoffDistance)
                        + chargeDipole(ForceOrPotential::Potential, atomConfiguration, cutoffDistance);

                double potentialShifted = potential - potentialAtRCut + (r - cutoffDistance)*forceAtRCut;
                double forceShifted = force - forceAtRCut;

                // Put this value into precomputed table
                at(at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1),mappedAtomicNumber2),i) = forceShifted;
                at(at(at(m_precomputedTwoParticlePotential,mappedAtomicNumber1),mappedAtomicNumber2),i) = potentialShifted;
            }
        }

    }
}

double USCWaterPotentialParameters::stericRepulsion(ForceOrPotential forceOrPotential, int atomConfiguration, double r) {
    if(forceOrPotential == ForceOrPotential::Force) return at(H_ij,atomConfiguration)*at(eta_ij,atomConfiguration)*pow(r, -at(eta_ij,atomConfiguration) - 1);
    else return at(H_ij,atomConfiguration)*pow(r, -at(eta_ij, atomConfiguration));
}

double USCWaterPotentialParameters::coulomb(ForceOrPotential forceOrPotential, int atomConfiguration, int atomicNumber1, int atomicNumber2, double r) {
    if(forceOrPotential == ForceOrPotential::Force) return at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(r, -2)*exp(-r/r1s)
            + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(r, -1)*exp(-r/r1s)/r1s;

    else return at(Z_i, atomicNumber1)*at(Z_i, atomicNumber2)*pow(r, -1)*exp(-r/r1s);
}

double USCWaterPotentialParameters::chargeDipole(ForceOrPotential forceOrPotential, int atomConfiguration, double r) {
    // Wrong units in the paper, using these values instead.
    if(forceOrPotential == ForceOrPotential::Force) return -at(D_ij,atomConfiguration)*0.5*4*pow(r, -5)*exp(-r/at(r4s,atomConfiguration)) - at(D_ij,atomConfiguration)*0.5*pow(r, -4)*exp(-r/at(r4s,atomConfiguration))/at(r4s,atomConfiguration);
    else return -0.5*at(D_ij,atomConfiguration)*pow(r, -4)*exp(-r/at(r4s, atomConfiguration));
}
