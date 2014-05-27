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
    atomManager.setCutoffDistance(2*m_maxTwoParticleCutoffDistance);
//    m_iteratorDefault.setMaximumNeighborDistance(2*m_maxThreeParticleCutoffDistance);
//    m_iteratorDefault.iterate(atomManager);
    m_iteratorAllPairs.iterate(atomManager);
}

// #define PRECOMPUTED
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

#ifdef PRECOMPUTED
    int mappedAtomicNumber1 = at(m_atomicNumberMap, atomicNumber1);
    int mappedAtomicNumber2 = at(m_atomicNumberMap, atomicNumber2);
    int precomputedTableIndex = r2*at(oneOverCutoffDistancesSquared, atomConfiguration)*m_numberOfPrecomputedTwoParticleForces;

    // This test is the same as if(r2 > cutoffDistnaceSquared)
    if(precomputedTableIndex>=m_numberOfPrecomputedTwoParticleForces) return;

    // Forces in this and the next bin to interpolate
    vector<double> &table = at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1),mappedAtomicNumber2);
    double force0 = at(table,precomputedTableIndex);
    double force1 = at(table,precomputedTableIndex+1);

    // Linearly interpolate between these values
    double force = force0 + (force1 - force0)*(r2 - precomputedTableIndex*m_deltaR2)*m_oneOverDeltaR2;

    double r = sqrt(r2);
    double oneOverR = 1.0/r;
    double oneOverR2 = oneOverR*oneOverR;
    double oneOverR3 = oneOverR2*oneOverR;
    double oneOverR4 = oneOverR2*oneOverR2;

#else

    if(r2 > at(cutoffDistancesSquared,atomConfiguration)) return;
#ifdef DEBUG
    if(r2 < 1e-5) {
        std::cout << "Error in USCSIO2Potential::twoParticleAction. Relative distance is near zero:" << std::endl;
        std::cout << *atom1 << endl;
        std::cout << *atom2 << endl;
        exit(1);
    }
#endif
    double r = sqrt(r2);
    double oneOverR = 1.0/r;
    double oneOverR2 = oneOverR*oneOverR;
    double oneOverR3 = oneOverR2*oneOverR;
    double oneOverR4 = oneOverR2*oneOverR2;
    double oneOverR5 = oneOverR2*oneOverR3;
    double oneOverR6 = oneOverR3*oneOverR3;

    double force = at(H_ij,atomConfiguration)*at(eta_ij,atomConfiguration)*pow(r, -at(eta_ij,atomConfiguration) - 2)
                    + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*oneOverR3*exp(-r*at(oneOverR1s,atomConfiguration))
                    + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*oneOverR2*exp(-r*at(oneOverR1s,atomConfiguration))*at(oneOverR1s,atomConfiguration)
                    - at(D_ij,atomConfiguration)*0.5*4*oneOverR6*exp(-r*at(oneOverR4s,atomConfiguration))
                    - at(D_ij,atomConfiguration)*0.5*oneOverR5*exp(-r*at(oneOverR4s,atomConfiguration))*at(oneOverR4s,atomConfiguration);
#endif

    double potentialEnergy = at(H_ij,atomConfiguration)*pow(oneOverR,at(eta_ij,atomConfiguration))
                             + Z_i.at(atomicNumber1)*Z_i.at(atomicNumber2)*oneOverR*exp(-r*at(oneOverR1s,atomConfiguration))
                             - 0.5*at(D_ij,atomConfiguration)*oneOverR4*exp(-r*at(oneOverR4s,atomConfiguration));
    int numberOfGhosts = atom1->ghost() + atom2->ghost();

    m_potentialEnergy += 0.5*potentialEnergy*(2-numberOfGhosts);


    atom1->force[0] += force*dx;
    atom1->force[1] += force*dy;
    atom1->force[2] += force*dz;

    atom2->force[0] -= force*dx;
    atom2->force[1] -= force*dy;
    atom2->force[2] -= force*dz;
}

void USCSIO2Potential::threeParticleAction(Atom *atomi, Atom *atomj, Atom *atomk)
{
    int atomicNumber1 = atomi->type()->atomicNumber();
    int atomicNumber2 = atomj->type()->atomicNumber();
    int atomicNumber3 = atomk->type()->atomicNumber();
    int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),atomicNumber3);
    if(atomConfiguration == +AtomConfiguration::NotUsed) return; // This configuration has zero contribution to the force
    sortAtoms(atomj, atomi, atomk, atomConfiguration);
    atomicNumber1 = atomi->type()->atomicNumber();
    atomicNumber2 = atomj->type()->atomicNumber();
    atomicNumber3 = atomk->type()->atomicNumber();

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
    double oneOverCosThetaIJKMinusCosThetaZero = 1.0/cosThetaIJKMinusCosThetaZero;

    double Vijk = at(B_ijk,atomConfiguration)
                  *exp(at(ksi,atomConfiguration)*(oneOverRijMinusRzero + oneOverRikMinusRzero))
                  *pow(cosThetaIJKMinusCosThetaZero, 2);

    int numberOfGhosts = atomi->ghost() + atomj->ghost() + atomk->ghost();
    m_potentialEnergy += 0.3333333333*Vijk*(3-numberOfGhosts);

    // Potential derivatives
    double dVijk_dCosThetaijk = Vijk*2*oneOverCosThetaIJKMinusCosThetaZero;
    double dVijk_dRij = -Vijk*at(ksi,atomConfiguration)*pow(oneOverRijMinusRzero, 2);
    double dVijk_dRik = -Vijk*at(ksi,atomConfiguration)*pow(oneOverRikMinusRzero, 2);

    // Theta derivatives
    // X
    double dCosThetaijk_dXij = oneOverRij*oneOverRik*(xik + xij*rijDotRik*oneOverRij*oneOverRij);
    double dCosThetaijk_dXik = oneOverRij*oneOverRik*(xij + xik*rijDotRik*oneOverRik*oneOverRik);
    // Y
    double dCosThetaijk_dYij = oneOverRij*oneOverRik*(yik + yij*rijDotRik*oneOverRij*oneOverRij);
    double dCosThetaijk_dYik = oneOverRij*oneOverRik*(yij + yik*rijDotRik*oneOverRik*oneOverRik);
    // Z
    double dCosThetaijk_dZij = oneOverRij*oneOverRik*(zik + zij*rijDotRik*oneOverRij*oneOverRij);
    double dCosThetaijk_dZik = oneOverRij*oneOverRik*(zij + zik*rijDotRik*oneOverRik*oneOverRik);

    // RijRikDerivatives
    // X
    double dRij_dXij = xij*oneOverRij;
    double dRik_dXik = xik*oneOverRik;
    // Y
    double dRij_dYij = yij*oneOverRij;
    double dRik_dYik = yik*oneOverRik;
    // Z
    double dRij_dZij = zij*oneOverRij;
    double dRik_dZik = zik*oneOverRik;

    // XijXikDerivatives
    // X
    double dXij_xi = 1;
    double dXik_xi = 1;
    double dXij_xj = -1;
    double dXik_xk = -1;
    double dXik_xj = 0;
    double dXij_xk = 0;

    // Y
    double dYij_yi = 1;
    double dYik_yi = 1;
    double dYij_yj = -1;
    double dYik_yk = -1;
    double dYik_yj = 0;
    double dYij_yk = 0;
    // Z
    double dZij_zi = 1;
    double dZik_zi = 1;
    double dZij_zj = -1;
    double dZik_zk = -1;
    double dZik_zj = 0;
    double dZij_zk = 0;

    double Fxi = -(dVijk_dRij*dRij_dXij*dXij_xi + dVijk_dRik*dRik_dXik*dXik_xi + dVijk_dCosThetaijk*dCosThetaijk_dXij*dXij_xi + dVijk_dCosThetaijk*dCosThetaijk_dXik*dXik_xi);
    double Fyi = -(dVijk_dRij*dRij_dYij*dYij_yi + dVijk_dRik*dRik_dYik*dYik_yi + dVijk_dCosThetaijk*dCosThetaijk_dYij*dYij_yi + dVijk_dCosThetaijk*dCosThetaijk_dYik*dYik_yi);
    double Fzi = -(dVijk_dRij*dRij_dZij*dZij_zi + dVijk_dRik*dRik_dZik*dZik_zi + dVijk_dCosThetaijk*dCosThetaijk_dZij*dZij_zi + dVijk_dCosThetaijk*dCosThetaijk_dZik*dZik_zi);

    double Fxj = -(dVijk_dRij*dRij_dXij*dXij_xj + dVijk_dRik*dRik_dXik*dXik_xj + dVijk_dCosThetaijk*dCosThetaijk_dXij*dXij_xj + dVijk_dCosThetaijk*dCosThetaijk_dXik*dXik_xj);
    double Fyj = -(dVijk_dRij*dRij_dYij*dYij_yj + dVijk_dRik*dRik_dYik*dYik_yj + dVijk_dCosThetaijk*dCosThetaijk_dYij*dYij_yj + dVijk_dCosThetaijk*dCosThetaijk_dYik*dYik_yj);
    double Fzj = -(dVijk_dRij*dRij_dZij*dZij_zj + dVijk_dRik*dRik_dZik*dZik_zj + dVijk_dCosThetaijk*dCosThetaijk_dZij*dZij_zj + dVijk_dCosThetaijk*dCosThetaijk_dZik*dZik_zj);

    double Fxk = -(dVijk_dRij*dRij_dXij*dXij_xk + dVijk_dRik*dRik_dXik*dXik_xk + dVijk_dCosThetaijk*dCosThetaijk_dXij*dXij_xk + dVijk_dCosThetaijk*dCosThetaijk_dXik*dXik_xk);
    double Fyk = -(dVijk_dRij*dRij_dYij*dYij_yk + dVijk_dRik*dRik_dYik*dYik_yk + dVijk_dCosThetaijk*dCosThetaijk_dYij*dYij_yk + dVijk_dCosThetaijk*dCosThetaijk_dYik*dYik_yk);
    double Fzk = -(dVijk_dRij*dRij_dZij*dZij_zk + dVijk_dRik*dRik_dZik*dZik_zk + dVijk_dCosThetaijk*dCosThetaijk_dZij*dZij_zk + dVijk_dCosThetaijk*dCosThetaijk_dZik*dZik_zk);

    atomi->force[0] += Fxi;
    atomi->force[1] += Fyi;
    atomi->force[2] += Fzi;

    atomj->force[0] += Fxj;
    atomj->force[1] += Fyj;
    atomj->force[2] += Fzj;

    atomk->force[0] += Fxk;
    atomk->force[1] += Fyk;
    atomk->force[2] += Fzk;
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

    D_ij.at(+AtomConfiguration::Si_O) = 3.456*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),4);
    D_ij.at(+AtomConfiguration::Si_Si) = 0.0;
    D_ij.at(+AtomConfiguration::O_O) = 1.728*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),4);

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
    double rMinSquared = 0;
    double rMaxSquared = m_maxTwoParticleCutoffDistance*m_maxTwoParticleCutoffDistance;
    m_deltaR2 = (rMaxSquared - rMinSquared) / (m_numberOfPrecomputedTwoParticleForces-1);
    m_oneOverDeltaR2 = 1.0/m_deltaR2;

    for(int atomicNumber1 : m_activeAtomTypes) {
        int mappedAtomicNumber1 = m_atomicNumberMap.at(atomicNumber1);
        for(int atomicNumber2 : m_activeAtomTypes) {
            int mappedAtomicNumber2 = m_atomicNumberMap.at(atomicNumber2);
            for(int i=0; i<m_numberOfPrecomputedTwoParticleForces; i++) {
                double r2 = rMinSquared + i*m_deltaR2;

                int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),+AtomTypes::NoAtom);
                if(atomConfiguration == +AtomConfiguration::NotUsed) continue; // This configuration has zero contribution to the force
                if(r2 > at(cutoffDistancesSquared,atomConfiguration)) continue;

                double r = sqrt(r2);
                double oneOverR = 1.0/r;
                double oneOverR2 = oneOverR*oneOverR;
                double oneOverR3 = oneOverR2*oneOverR;
                double oneOverR5 = oneOverR2*oneOverR3;
                double oneOverR6 = oneOverR3*oneOverR3;

                double force = at(H_ij,atomConfiguration)*at(eta_ij,atomConfiguration)*pow(r, -at(eta_ij,atomConfiguration) - 2)
                                + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*oneOverR3*exp(-r*at(oneOverR1s,atomConfiguration))
                                + at(Z_i,atomicNumber1)*at(Z_i,atomicNumber2)*oneOverR2*exp(-r*at(oneOverR1s,atomConfiguration))*at(oneOverR1s,atomConfiguration)
                                - at(D_ij,atomConfiguration)*0.5*4*oneOverR6*exp(-r*at(oneOverR4s,atomConfiguration))
                                - at(D_ij,atomConfiguration)*0.5*oneOverR5*exp(-r*at(oneOverR4s,atomConfiguration))*at(oneOverR4s,atomConfiguration);

                // Put this value into precomputed table
                at(at(at(m_precomputedTwoParticleForces,mappedAtomicNumber1),mappedAtomicNumber2),i) = force;
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
