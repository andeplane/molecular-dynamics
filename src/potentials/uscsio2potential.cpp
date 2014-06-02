/* USC SiO2 potential
 * [1] Vashishta, P., et al. "Interaction potential for SiO 2: a molecular-dynamics study of structural correlations." Physical Review B 41.17 (1990): 12197.
 */

#include <potentials/uscsio2potential.h>
#include <unitconverter.h>
#include <cmath>
#include <includes.h>
using CompPhys::Utils::at;

void USCSIO2Potential::calculateForces(AtomManager &atomManager)
{
    m_potentialEnergy = 0;
    double maxCutoffDistance = 2*m_maxThreeParticleCutoffDistance;
    maxCutoffDistance = std::max(maxCutoffDistance, m_maxTwoParticleCutoffDistance);

    atomManager.setCutoffDistance(maxCutoffDistance);
    m_iteratorDefault.setMaximumNeighborDistance(2*m_maxThreeParticleCutoffDistance);
    m_iteratorDefault.iterate(atomManager);
}

void USCSIO2Potential::twoParticleAction(Atom *atom1, Atom *atom2)
{
    int atomicNumber1 = atom1->type()->atomicNumber();
    int atomicNumber2 = atom2->type()->atomicNumber();
    int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),+AtomTypes::NoAtom);

    if(atomConfiguration == +AtomConfiguration::NotUsed) return; // This configuration has zero contribution to the force

    double dx = atom1->position[0] - atom2->position[0];
    double dy = atom1->position[1] - atom2->position[1];
    double dz = atom1->position[2] - atom2->position[2];
    double r2 = dx*dx + dy*dy + dz*dz;

    int mappedAtomicNumber1 = at(m_atomicNumberMap, atomicNumber1);
    int mappedAtomicNumber2 = at(m_atomicNumberMap, atomicNumber2);

    int precomputedTableIndex = r2*at(oneOverCutoffDistancesSquared, atomConfiguration)*m_numberOfPrecomputedTwoParticleForces;

    // This test is the same as if(r2 > cutoffDistanceSquared)
    if(precomputedTableIndex>=m_numberOfPrecomputedTwoParticleForces) return;

    // Forces in this and the next bin to interpolate
    vector<double> &forceTable = at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1),mappedAtomicNumber2);
    vector<double> &energyTable = at(at(m_precomputedTwoParticlePotential,mappedAtomicNumber1),mappedAtomicNumber2);
    double force0 = at(forceTable,precomputedTableIndex);
    double force1 = at(forceTable,precomputedTableIndex+1);

    double energy0 = at(energyTable,precomputedTableIndex);
    double energy1 = at(energyTable,precomputedTableIndex+1);

    // Linearly interpolate between these values
    double force = force0 + (force1 - force0)*(r2 - precomputedTableIndex*m_deltaR2)*m_oneOverDeltaR2;
    double potentialEnergy = energy0 + (energy1 - energy0)*(r2 - precomputedTableIndex*m_deltaR2)*m_oneOverDeltaR2;

    int numberOfGhosts = atom1->ghost() + atom2->ghost();

    m_potentialEnergy += 0.5*potentialEnergy*(2-numberOfGhosts);

    double oneOverR = 1.0/sqrt(r2);
    atom1->force[0] += force*dx*oneOverR;
    atom1->force[1] += force*dy*oneOverR;
    atom1->force[2] += force*dz*oneOverR;

    atom2->force[0] -= force*dx*oneOverR;
    atom2->force[1] -= force*dy*oneOverR;
    atom2->force[2] -= force*dz*oneOverR;
}

double dVdXij(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*xij*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2)), 2)*sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(-2*xij*(xij*xik + yij*yik + zij*zik)/(pow(pow(xij, 2) + pow(yij, 2) + pow(zij, 2), 3.0L/2.0L)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + 2*xik/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

double dVdYij(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*yij*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2)), 2)*sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(-2*yij*(xij*xik + yij*yik + zij*zik)/(pow(pow(xij, 2) + pow(yij, 2) + pow(zij, 2), 3.0L/2.0L)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + 2*yik/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

double dVdZij(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*zij*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2)), 2)*sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(-2*zij*(xij*xik + yij*yik + zij*zik)/(pow(pow(xij, 2) + pow(yij, 2) + pow(zij, 2), 3.0L/2.0L)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + 2*zik/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

double dVdXik(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*xik*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2)), 2)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(2*xij/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) - 2*xik*(xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*pow(pow(xik, 2) + pow(yik, 2) + pow(zik, 2), 3.0L/2.0L)))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

double dVdYik(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*yik*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2)), 2)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(2*yij/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) - 2*yik*(xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*pow(pow(xik, 2) + pow(yik, 2) + pow(zik, 2), 3.0L/2.0L)))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

double dVdZik(double B, double costheta0, double ksi, double r0, double xij, double xik, double yij, double yik, double zij, double zik) {
   return -B*ksi*zik*pow(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))), 2)*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))))/(pow(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2)), 2)*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + B*(-costheta0 + (xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))))*(2*zij/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) - 2*zik*(xij*xik + yij*yik + zij*zik)/(sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))*pow(pow(xik, 2) + pow(yik, 2) + pow(zik, 2), 3.0L/2.0L)))*exp(ksi/(-r0 + sqrt(pow(xik, 2) + pow(yik, 2) + pow(zik, 2))) + ksi/(-r0 + sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2))));
}

void USCSIO2Potential::threeParticleAction(Atom *atomi, Atom *atomj, Atom *atomk)
{
    // return;
    int atomicNumber1 = atomi->type()->atomicNumber();
    int atomicNumber2 = atomj->type()->atomicNumber();
    int atomicNumber3 = atomk->type()->atomicNumber();
    int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),atomicNumber3);
    if(atomConfiguration == +AtomConfiguration::NotUsed) return; // This configuration has zero contribution to the force

    sortAtoms(atomj, atomi, atomk, atomConfiguration);

    double xij = atomi->position[0] - atomj->position[0];
    double xik = atomi->position[0] - atomk->position[0];

    double yij = atomi->position[1] - atomj->position[1];
    double yik = atomi->position[1] - atomk->position[1];

    double zij = atomi->position[2] - atomj->position[2];
    double zik = atomi->position[2] - atomk->position[2];

    double rij = sqrt(xij*xij + yij*yij + zij*zij);
    double rik = sqrt(xik*xik + yik*yik + zik*zik);

    if(rij > at(r0,atomConfiguration) || rik > at(r0,atomConfiguration)) return;

    double rijDotRik = xij*xik + yij*yik + zij*zik;

    double oneOverRij = 1.0/rij;
    double oneOverRik = 1.0/rik;
    double oneOverRijMinusRzero = 1.0/(rij - at(r0,atomConfiguration));
    double oneOverRikMinusRzero = 1.0/(rik - at(r0,atomConfiguration));
    double cosThetaIJK = rijDotRik*oneOverRij*oneOverRik;
    double cosThetaIJKMinusCosThetaZero = cosThetaIJK - at(cosThetaZero,atomConfiguration);

    double Vijk = at(B_ijk,atomConfiguration)
                  *exp(at(ksi,atomConfiguration)*(oneOverRijMinusRzero + oneOverRikMinusRzero))
                  *pow(cosThetaIJKMinusCosThetaZero, 2);

    int numberOfGhosts = atomi->ghost() + atomj->ghost() + atomk->ghost();
    m_potentialEnergy += 0.3333333333*Vijk*(3-numberOfGhosts);

    double B_ijk = at(this->B_ijk,atomConfiguration);
    double cosThetaZero = at(this->cosThetaZero, atomConfiguration);
    double ksi = at(this->ksi, atomConfiguration);
    double r0 = at(this->r0, atomConfiguration);
    double dVdXij2 = dVdXij(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);
    double dVdYij2 = dVdYij(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);
    double dVdZij2 = dVdZij(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);

    double dVdXik2 = dVdXik(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);
    double dVdYik2 = dVdYik(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);
    double dVdZik2 = dVdZik(B_ijk, cosThetaZero, ksi, r0, xij, xik, yij, yik, zij, zik);

    atomi->force[0] -= dVdXij2 + dVdXik2;
    atomi->force[1] -= dVdYij2 + dVdYik2;
    atomi->force[2] -= dVdZij2 + dVdZik2;

    atomj->force[0] += dVdXij2;
    atomj->force[1] += dVdYij2;
    atomj->force[2] += dVdZij2;

    atomk->force[0] += dVdXik2;
    atomk->force[1] += dVdYik2;
    atomk->force[2] += dVdZik2;
}

std::string USCSIO2Potential::coefficientString() const
{
    char output[10000];
    sprintf(output, "                     Si                     O                          \n");
    sprintf(output,"%s Zi [e]             %.2f                  %.2f                         \n",output,Z_i.at(+AtomTypes::Silicon),Z_i.at(+AtomTypes::Oxygen));
    sprintf(output,"%sr1s=%.2fÅ         r4s=%.2fÅ              rc=%.2fÅ       e=%eC         \n",output,UnitConverter::lengthToAngstroms(r1s.at(+AtomConfiguration::Si_O)), UnitConverter::lengthToAngstroms(r4s.at(+AtomConfiguration::Si_O)), UnitConverter::lengthToAngstroms(cutoffDistances.at(+AtomConfiguration::Si_O)), UnitConverter::chargeToSI(1.0));
    sprintf(output,"%s                              Two body                                 \n",output);
    sprintf(output,"%s                   Si-Si                  Si-O                  O-O    \n",output);
    sprintf(output,"%seta_ij              %d                     %d                     %d    \n",output, eta_ij.at(+AtomConfiguration::Si_Si), eta_ij.at(+AtomConfiguration::Si_O), eta_ij.at(+AtomConfiguration::O_O));
    sprintf(output,"%sH_ij [eVÅ^eta_ij] %f              %f            %f    \n", output,H_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_Si)), H_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_O)), H_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O)));
    sprintf(output,"%sD_ij [eVÅ^4]      %f              %f               %f    \n", output,D_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4), D_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4), D_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4));
    sprintf(output,"%sW_ij [eVÅ^6]      %f              %f               %f    \n", output,W_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6), W_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6), W_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6));
    sprintf(output,"%s                             Three body                                 \n",output);
    sprintf(output,"%s         Bijk [eV]    theta_0 [deg]    Cijk    ksi [Å]    r0 [Å]          \n",output);
    sprintf(output,"%sO-Si-O   %.3f        %.2f           %.2f    %.2f       %.2f               \n", output, UnitConverter::energyToEv(B_ijk.at(+AtomConfiguration::O_Si_O)), UnitConverter::radiansToDegrees(acos(cosThetaZero.at(+AtomConfiguration::O_Si_O))), C_ijk.at(+AtomConfiguration::O_Si_O), UnitConverter::lengthToAngstroms(ksi.at(+AtomConfiguration::O_Si_O)), UnitConverter::lengthToAngstroms(r0.at(+AtomConfiguration::O_Si_O)));
    sprintf(output,"%sSi-O-Si  %.3f       %.2f           %.2f    %.2f       %.2f               \n", output,UnitConverter::energyToEv(B_ijk.at(+AtomConfiguration::Si_O_Si)), UnitConverter::radiansToDegrees(acos(cosThetaZero.at(+AtomConfiguration::Si_O_Si))), C_ijk.at(+AtomConfiguration::Si_O_Si), UnitConverter::lengthToAngstroms(ksi.at(+AtomConfiguration::Si_O_Si)), UnitConverter::lengthToAngstroms(r0.at(+AtomConfiguration::Si_O_Si)));

    return string(output);
}


int USCSIO2Potential::numberOfPrecomputedTwoParticleForces() const
{
    return m_numberOfPrecomputedTwoParticleForces;
}

void USCSIO2Potential::setNumberOfPrecomputedTwoParticleForces(int numberOfPrecomputedTwoParticleForces)
{
    m_numberOfPrecomputedTwoParticleForces = numberOfPrecomputedTwoParticleForces;
    calculatePrecomputedTwoParticleForces();
}

void USCSIO2Potential::initialize() {
    int highestAtomicNumber = +AtomTypes::Silicon;
    int numberOfAtomicNumbers = highestAtomicNumber+1;
    UnitConverter uc;

    m_activeAtomTypes.push_back(+AtomTypes::Oxygen);
    m_activeAtomTypes.push_back(+AtomTypes::Silicon);
    m_atomicNumberMap.resize(numberOfAtomicNumbers,-1);
    int atomMapIndex = 0; // First used atom will be the first index
    for(int atomType : m_activeAtomTypes) {
        m_atomicNumberMap.at(atomType) = atomMapIndex++;
    }

    B_ijk.resize(+AtomConfiguration::NumberOfConfigurations,0);
    cosThetaZero.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r0.resize(+AtomConfiguration::NumberOfConfigurations,0);
    ksi.resize(+AtomConfiguration::NumberOfConfigurations,0);
    C_ijk.resize(+AtomConfiguration::NumberOfConfigurations,0);
    eta_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    H_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    D_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    W_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    cutoffDistances.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r1s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r4s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    oneOverR1s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    oneOverR4s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    Z_i.resize(numberOfAtomicNumbers,0);
    alpha_i.resize(numberOfAtomicNumbers,0);

    // ALL THESE VALUES ARE FROM TABLE 1 IN [1]
    r1s.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(4.43);
    r1s.at(+AtomConfiguration::Si_Si) = uc.lengthFromAngstroms(4.43);
    r1s.at(+AtomConfiguration::O_O) = uc.lengthFromAngstroms(4.43);

    r4s.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::Si_Si) = uc.lengthFromAngstroms(2.5);
    r4s.at(+AtomConfiguration::O_O) = uc.lengthFromAngstroms(2.5);

    oneOverR1s.at(+AtomConfiguration::Si_O) = 1.0/r1s.at(+AtomConfiguration::Si_O);
    oneOverR1s.at(+AtomConfiguration::Si_Si) = 1.0/r1s.at(+AtomConfiguration::Si_Si);
    oneOverR1s.at(+AtomConfiguration::O_O) = 1.0/r1s.at(+AtomConfiguration::O_O);

    oneOverR4s.at(+AtomConfiguration::Si_O) = 1.0/r4s.at(+AtomConfiguration::Si_O);
    oneOverR4s.at(+AtomConfiguration::Si_Si) = 1.0/r4s.at(+AtomConfiguration::Si_Si);
    oneOverR4s.at(+AtomConfiguration::O_O) = 1.0/r4s.at(+AtomConfiguration::O_O);

    Z_i.at(+AtomTypes::Silicon) = 1.2; // Given in atomic units
    Z_i.at(+AtomTypes::Oxygen) = -0.6; // Given in atomic units
    alpha_i.at(+AtomTypes::Silicon) = 0;
    alpha_i.at(+AtomTypes::Oxygen) = 2.4*pow(UnitConverter::lengthFromAngstroms(1.0),3);

    eta_ij.at(+AtomConfiguration::Si_O) = 9;
    eta_ij.at(+AtomConfiguration::Si_Si) = 11;
    eta_ij.at(+AtomConfiguration::O_O) = 7;

    cutoffDistances.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(5.5);
    cutoffDistances.at(+AtomConfiguration::Si_Si) = uc.lengthFromAngstroms(5.5);
    cutoffDistances.at(+AtomConfiguration::O_O) = uc.lengthFromAngstroms(5.5);

    // Copy the values and square them
    cutoffDistancesSquared = cutoffDistances;
    for(double &cutoffDistanceSquaredValue : cutoffDistancesSquared) {
        cutoffDistanceSquaredValue = cutoffDistanceSquaredValue*cutoffDistanceSquaredValue;
    }

    // Copy the values and take inverse
    oneOverCutoffDistancesSquared = cutoffDistancesSquared;
    for(double &oneOverCutoffDistanceSquaredValue : oneOverCutoffDistancesSquared) {
        oneOverCutoffDistanceSquaredValue = 1.0/oneOverCutoffDistanceSquaredValue;
    }

    D_ij.at(+AtomConfiguration::Si_O) = 3.456*pow(uc.lengthFromAngstroms(1.0),3);
    D_ij.at(+AtomConfiguration::Si_Si) = 0.0;
    D_ij.at(+AtomConfiguration::O_O) = 1.728*   pow(uc.lengthFromAngstroms(1.0),3);

    H_ij.at(+AtomConfiguration::Si_O) = 78.3143*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_O));
    H_ij.at(+AtomConfiguration::Si_Si) = 0.39246*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_Si));
    H_ij.at(+AtomConfiguration::O_O) = 355.5263*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O));

    B_ijk.at(+AtomConfiguration::O_Si_O) = uc.energyFromEv(4.993);
    B_ijk.at(+AtomConfiguration::Si_O_Si) = uc.energyFromEv(19.972);

    cosThetaZero.at(+AtomConfiguration::O_Si_O) = cos(uc.degreesToRadians(109.47));
    cosThetaZero.at(+AtomConfiguration::Si_O_Si) = cos(uc.degreesToRadians(141.0));

    r0.at(+AtomConfiguration::O_Si_O) = uc.lengthFromAngstroms(2.60);
    r0.at(+AtomConfiguration::Si_O_Si) = uc.lengthFromAngstroms(2.60);

    ksi.at(+AtomConfiguration::O_Si_O) = uc.lengthFromAngstroms(1.0);
    ksi.at(+AtomConfiguration::Si_O_Si) = uc.lengthFromAngstroms(1.0);

    C_ijk.at(+AtomConfiguration::O_Si_O) = 0;
    C_ijk.at(+AtomConfiguration::Si_O_Si) = 0;

    // Remember the maximum cutoff distance. This will be used to generate the cells before two particle forces
    m_maxTwoParticleCutoffDistance = *std::max_element(cutoffDistances.begin(), cutoffDistances.end());
    m_maxThreeParticleCutoffDistance = *std::max_element(r0.begin(), r0.end());

    m_configurationMap.resize(numberOfAtomicNumbers,vector<vector<int> >(numberOfAtomicNumbers, vector<int>(numberOfAtomicNumbers,0)));
    // Two particle configurations
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen).at(+AtomTypes::NoAtom) = +AtomConfiguration::O_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_O;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Silicon).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_Si;

    // Si-O-Si
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_Si;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_Si;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen) = +AtomConfiguration::Si_O_Si;

    // O-Si-O
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon) = +AtomConfiguration::O_Si_O;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen) = +AtomConfiguration::O_Si_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen) = +AtomConfiguration::O_Si_O;
}

void USCSIO2Potential::calculatePrecomputedTwoParticleForces()
{
    m_precomputedTwoParticleForces.resize(m_activeAtomTypes.size(),
                                          vector<vector<double > >(m_activeAtomTypes.size(),
                                          vector<double >(m_numberOfPrecomputedTwoParticleForces+1,0))); // +1 because we need rom for the extra value to interpolate in the last bin

    m_precomputedTwoParticlePotential.resize(m_activeAtomTypes.size(),
                                          vector<vector<double > >(m_activeAtomTypes.size(),
                                          vector<double >(m_numberOfPrecomputedTwoParticleForces+1,0))); // +1 because we need rom for the extra value to interpolate in the last bin


    double rMinSquared = 0;
    double rMaxSquared = m_maxTwoParticleCutoffDistance*m_maxTwoParticleCutoffDistance;
    m_deltaR2 = (rMaxSquared - rMinSquared) / (m_numberOfPrecomputedTwoParticleForces-1);
    m_oneOverDeltaR2 = 1.0/m_deltaR2;

    for(int atomicNumber1 : m_activeAtomTypes) {
        int mappedAtomicNumber1 = m_atomicNumberMap.at(atomicNumber1);
        for(int atomicNumber2 : m_activeAtomTypes) {
            int mappedAtomicNumber2 = m_atomicNumberMap.at(atomicNumber2);
            int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),+AtomTypes::NoAtom);
            if(atomConfiguration == +AtomConfiguration::NotUsed) continue; // This configuration has zero contribution to the force

            double Dij = at(alpha_i,atomicNumber1)*at(Z_i,atomicNumber2)*at(Z_i,atomicNumber2) + at(alpha_i,atomicNumber2)*at(Z_i,atomicNumber1)*at(Z_i,atomicNumber1);
//            cout << "D_ij (new paper: " << at(D_ij, atomConfiguration) << endl;
//            cout << "Dij (old paper: " << Dij << endl;
//            cout << "Ratio: " << Dij / at(D_ij, atomConfiguration) << endl;

            for(int i=0; i<m_numberOfPrecomputedTwoParticleForces; i++) {
                double r2 = rMinSquared + i*m_deltaR2;

                if(r2 > at(cutoffDistancesSquared,atomConfiguration)) continue;
                double rCut = at(cutoffDistances, atomConfiguration);
                double r = sqrt(r2);

                double force = at(H_ij,atomConfiguration)*at(eta_ij,atomConfiguration)*pow(r, -at(eta_ij,atomConfiguration) - 1)
                                + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(r, -2)*exp(-r*at(oneOverR1s,atomConfiguration))
                                + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(r, -1)*exp(-r*at(oneOverR1s,atomConfiguration))*at(oneOverR1s,atomConfiguration)
                                - Dij*0.5*4*pow(r, -5)*exp(-r*at(oneOverR4s,atomConfiguration))
                                - Dij*0.5*pow(r, -4)*exp(-r*at(oneOverR4s,atomConfiguration))*at(oneOverR4s,atomConfiguration);

                double forceAtRCut = at(H_ij,atomConfiguration)*at(eta_ij,atomConfiguration)*pow(rCut, -at(eta_ij,atomConfiguration) - 1)
                        + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(rCut, -2)*exp(-rCut*at(oneOverR1s,atomConfiguration))
                        + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*pow(rCut, -1)*exp(-rCut*at(oneOverR1s,atomConfiguration))*at(oneOverR1s,atomConfiguration)
                        - 0.5*Dij*4*pow(rCut, -5)*exp(-rCut*at(oneOverR4s,atomConfiguration))
                        - 0.5*Dij*pow(rCut, -4)*exp(-rCut*at(oneOverR4s,atomConfiguration))*at(oneOverR4s,atomConfiguration);

                double potential = at(H_ij,atomConfiguration)*pow(r, -at(eta_ij, atomConfiguration))
                                    + at(Z_i, atomicNumber1)*at(Z_i, atomicNumber2)*pow(r, -1)*exp(-r*at(oneOverR1s,atomConfiguration))
                                    - 0.5*Dij*pow(r, -4)*exp(-r*at(oneOverR4s, atomConfiguration));

                double potentialAtRCut = at(H_ij,atomConfiguration)*pow(rCut, -at(eta_ij, atomConfiguration))
                                    + at(Z_i, atomicNumber1)*at(Z_i, atomicNumber2)*pow(rCut, -1)*exp(-rCut*at(oneOverR1s,atomConfiguration))
                                    - 0.5*Dij*pow(rCut, -4)*exp(-rCut*at(oneOverR4s, atomConfiguration));

                double potentialShifted = potential - potentialAtRCut + (r - rCut)*forceAtRCut;
                double forceShifted = force - forceAtRCut;

                // Put this value into precomputed table
                at(at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1),mappedAtomicNumber2),i) = forceShifted;
                at(at(at(m_precomputedTwoParticlePotential,mappedAtomicNumber1),mappedAtomicNumber2),i) = potentialShifted;
            }
        }

    }

}

USCSIO2Potential::USCSIO2Potential()
{
    setName("USC SiO2");
    using namespace std::placeholders;
    m_iteratorDefault.setTwoParticleAction(std::bind(&USCSIO2Potential::twoParticleAction, this, _1, _2));
    m_iteratorDefault.setThreeParticleAction(std::bind(&USCSIO2Potential::threeParticleAction, this, _1, _2, _3));

    m_iteratorAllPairs.setTwoParticleAction(std::bind(&USCSIO2Potential::twoParticleAction, this, _1, _2));
    m_iteratorAllPairs.setThreeParticleAction(std::bind(&USCSIO2Potential::threeParticleAction, this, _1, _2, _3));

    initialize();
    setNumberOfPrecomputedTwoParticleForces(8192); // Default value
}

std::ostream& operator<<(std::ostream &stream, const USCSIO2Potential &potential) {
    return stream << potential.coefficientString().c_str();
}
