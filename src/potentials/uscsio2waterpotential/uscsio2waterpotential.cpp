/* USC SiO2 potential
 * [1] Vashishta, P., et al. "Interaction potential for SiO 2: a molecular-dynamics study of structural correlations." Physical Review B 41.17 (1990): 12197.
 */

#include <potentials/uscsio2waterpotential/uscsio2waterpotential.h>
#include <unitconverter.h>
#include <cmath>
#include <includes.h>
using CompPhys::Utils::at;

USCSIO2WaterPotential::USCSIO2WaterPotential()
{
    setName("USC SiO2-water");
    using namespace std::placeholders;
    m_iteratorDefault.setTwoParticleAction(std::bind(&USCSIO2WaterPotential::twoParticleAction, this, _1, _2));
    m_iteratorDefault.setThreeParticleAction(std::bind(&USCSIO2WaterPotential::threeParticleAction, this, _1, _2, _3));

    initialize();
    m_waterParameters = new USCWaterPotentialParameters(PotentialClass::Water, m_activeAtomTypes, m_atomicNumberMap, m_configurationMap);
    m_silicaParameters = new USCWaterPotentialParameters(PotentialClass::Silica, m_activeAtomTypes, m_atomicNumberMap, m_configurationMap);
    setNumberOfPrecomputedTwoParticleForces(8192); // Default value
}

void USCSIO2WaterPotential::calculateForces(AtomManager &atomManager)
{
    m_potentialEnergy = 0;
    double maxCutoffDistance = 2*m_waterParameters->cutoffDistance;

    atomManager.setCutoffDistance(m_twoParticleCutoffDistance);
    m_iteratorDefault.setMaximumNeighborDistance(2*m_threeParticleCutoffDistance);
    m_iteratorDefault.iterate(atomManager);
}

void USCSIO2WaterPotential::twoParticleAction(Atom *atom1, Atom *atom2)
{
    int atomicNumber1 = atom1->type()->atomicNumber();
    int atomicNumber2 = atom2->type()->atomicNumber();
    int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),+AtomTypes::NoAtom);

    if(atomConfiguration == +AtomConfiguration::NotUsed) return; // This configuration has zero contribution to the force

    double dx = atom1->position[0] - atom2->position[0];
    double dy = atom1->position[1] - atom2->position[1];
    double dz = atom1->position[2] - atom2->position[2];
    double r2 = dx*dx + dy*dy + dz*dz;

    int precomputedTableIndex = r2*m_twoParticleCutoffDistance*m_numberOfPrecomputedTwoParticleForces;

    // This test is the same as if(r2 > cutoffDistanceSquared)
    if(precomputedTableIndex>=m_numberOfPrecomputedTwoParticleForces) return;

    // Forces in this and the next bin to interpolate
    vector<double> &forceTable = m_silicaParameters->forcesTable(atomicNumber1, atomicNumber2);
    vector<double> &potentialTable = m_silicaParameters->potentialTable(atomicNumber1, atomicNumber2);

    double force0 = at(forceTable,precomputedTableIndex);
    double force1 = at(forceTable,precomputedTableIndex+1);

    double energy0 = at(potentialTable,precomputedTableIndex);
    double energy1 = at(potentialTable,precomputedTableIndex+1);

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

void USCSIO2WaterPotential::threeParticleAction(Atom *atomi, Atom *atomj, Atom *atomk)
{
    int atomicNumber1 = atomi->type()->atomicNumber();
    int atomicNumber2 = atomj->type()->atomicNumber();
    int atomicNumber3 = atomk->type()->atomicNumber();
    int atomConfiguration = at(at(at(m_configurationMap,atomicNumber1),atomicNumber2),atomicNumber3);
    if(atomConfiguration == +AtomConfiguration::NotUsed) return; // This configuration has zero contribution to the force

    // We know that the atoms are in ascending order (by atomId), now sort for correct configuration (Si-Si-O -> Si-O-Si)
    sortAtoms(atomj, atomi, atomk, atomConfiguration);

    double xij = atomi->position[0] - atomj->position[0];
    double xik = atomi->position[0] - atomk->position[0];

    double yij = atomi->position[1] - atomj->position[1];
    double yik = atomi->position[1] - atomk->position[1];

    double zij = atomi->position[2] - atomj->position[2];
    double zik = atomi->position[2] - atomk->position[2];

    double rij = sqrt(xij*xij + yij*yij + zij*zij);
    double rik = sqrt(xik*xik + yik*yik + zik*zik);

    double r0 = at(m_silicaParameters->r0, atomConfiguration);

    if(rij > r0 || rik > r0) return;

    double B = at(m_silicaParameters->B_ijk,atomConfiguration);
    double cosTheta0 = at(m_silicaParameters->cosThetaZero, atomConfiguration);
    double ksi = at(m_silicaParameters->ksi, atomConfiguration);

    double rijDotRik = xij*xik + yij*yik + zij*zik;
    double oneOverRij = 1.0/rij;
    double oneOverRik = 1.0/rik;
    double oneOverRijMinusRzero = 1.0/(rij - r0);
    double oneOverRikMinusRzero = 1.0/(rik - r0);
    double cosTheta = rijDotRik*oneOverRij*oneOverRik;
    double cosThetaMinusCosThetaZero = cosTheta - cosTheta0;

    double Vijk = B*exp(ksi*(oneOverRijMinusRzero + oneOverRikMinusRzero))
            *pow(cosThetaMinusCosThetaZero, 2);

    int numberOfGhosts = atomi->ghost() + atomj->ghost() + atomk->ghost();
    m_potentialEnergy += 0.3333333333*Vijk*(3-numberOfGhosts);

    // These expressions are from molecular-dynamics/python/threeParticleForce.py
    double dVdXij = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*xij*(cosTheta - cosTheta0) + 2*pow(r0 - rij, 2)*(cosTheta*rik*xij - rij*xik))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(pow(rij, 2)*rik*pow(r0 - rij, 2));
    double dVdYij = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*yij*(cosTheta - cosTheta0) + 2*pow(r0 - rij, 2)*(cosTheta*rik*yij - rij*yik))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(pow(rij, 2)*rik*pow(r0 - rij, 2));
    double dVdZij = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*zij*(cosTheta - cosTheta0) + 2*pow(r0 - rij, 2)*(cosTheta*rik*zij - rij*zik))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(pow(rij, 2)*rik*pow(r0 - rij, 2));

    double dVdXik = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*xik*(cosTheta - cosTheta0) + 2*pow(r0 - rik, 2)*(cosTheta*rij*xik - rik*xij))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(rij*pow(rik, 2)*pow(r0 - rik, 2));
    double dVdYik = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*yik*(cosTheta - cosTheta0) + 2*pow(r0 - rik, 2)*(cosTheta*rij*yik - rik*yij))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(rij*pow(rik, 2)*pow(r0 - rik, 2));
    double dVdZik = -B*(cosTheta - cosTheta0)*(ksi*rij*rik*zik*(cosTheta - cosTheta0) + 2*pow(r0 - rik, 2)*(cosTheta*rij*zik - rik*zij))*exp(ksi*(-2*r0 + rij + rik)/((r0 - rij)*(r0 - rik)))/(rij*pow(rik, 2)*pow(r0 - rik, 2));

    atomi->force[0] -= dVdXij + dVdXik;
    atomi->force[1] -= dVdYij + dVdYik;
    atomi->force[2] -= dVdZij + dVdZik;

    atomj->force[0] += dVdXij;
    atomj->force[1] += dVdYij;
    atomj->force[2] += dVdZij;

    atomk->force[0] += dVdXik;
    atomk->force[1] += dVdYik;
    atomk->force[2] += dVdZik;
}

std::string USCSIO2WaterPotential::coefficientString() const
{
    char output[10000];
//    sprintf(output, "                     Si                     O                          \n");
//    sprintf(output,"%s Zi [e]             %.2f                  %.2f                         \n",output,Z_i.at(+AtomTypes::Silicon),Z_i.at(+AtomTypes::Oxygen));
//    sprintf(output,"%sr1s=%.2fÅ         r4s=%.2fÅ              rc=%.2fÅ       e=%eC         \n",output,UnitConverter::lengthToAngstroms(r1s.at(+AtomConfiguration::Si_O)), UnitConverter::lengthToAngstroms(r4s.at(+AtomConfiguration::Si_O)), UnitConverter::lengthToAngstroms(cutoffDistances.at(+AtomConfiguration::Si_O)), UnitConverter::chargeToSI(1.0));
//    sprintf(output,"%s                              Two body                                 \n",output);
//    sprintf(output,"%s                   Si-Si                  Si-O                  O-O    \n",output);
//    sprintf(output,"%seta_ij              %d                     %d                     %d    \n",output, eta_ij.at(+AtomConfiguration::Si_Si), eta_ij.at(+AtomConfiguration::Si_O), eta_ij.at(+AtomConfiguration::O_O));
//    sprintf(output,"%sH_ij [eVÅ^eta_ij] %f              %f            %f    \n", output,H_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_Si)), H_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_O)), H_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O)));
//    sprintf(output,"%sD_ij [eVÅ^4]      %f              %f               %f    \n", output,D_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4), D_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4), D_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),4));
//    sprintf(output,"%sW_ij [eVÅ^6]      %f              %f               %f    \n", output,W_ij.at(+AtomConfiguration::Si_Si)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6), W_ij.at(+AtomConfiguration::Si_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6), W_ij.at(+AtomConfiguration::O_O)*UnitConverter::energyToEv(1.0)*pow(UnitConverter::lengthToAngstroms(1.0),6));
//    sprintf(output,"%s                             Three body                                 \n",output);
//    sprintf(output,"%s         Bijk [eV]    theta_0 [deg]    Cijk    ksi [Å]    r0 [Å]          \n",output);
//    sprintf(output,"%sO-Si-O   %.3f        %.2f           %.2f    %.2f       %.2f               \n", output, UnitConverter::energyToEv(B_ijk.at(+AtomConfiguration::O_Si_O)), UnitConverter::radiansToDegrees(acos(cosThetaZero.at(+AtomConfiguration::O_Si_O))), C_ijk.at(+AtomConfiguration::O_Si_O), UnitConverter::lengthToAngstroms(ksi.at(+AtomConfiguration::O_Si_O)), UnitConverter::lengthToAngstroms(r0.at(+AtomConfiguration::O_Si_O)));
//    sprintf(output,"%sSi-O-Si  %.3f       %.2f           %.2f    %.2f       %.2f               \n", output,UnitConverter::energyToEv(B_ijk.at(+AtomConfiguration::Si_O_Si)), UnitConverter::radiansToDegrees(acos(cosThetaZero.at(+AtomConfiguration::Si_O_Si))), C_ijk.at(+AtomConfiguration::Si_O_Si), UnitConverter::lengthToAngstroms(ksi.at(+AtomConfiguration::Si_O_Si)), UnitConverter::lengthToAngstroms(r0.at(+AtomConfiguration::Si_O_Si)));

    return string(output);
}


int USCSIO2WaterPotential::numberOfPrecomputedTwoParticleForces() const
{
    return m_numberOfPrecomputedTwoParticleForces;
}

void USCSIO2WaterPotential::setNumberOfPrecomputedTwoParticleForces(int numberOfPrecomputedTwoParticleForces)
{
    m_numberOfPrecomputedTwoParticleForces = numberOfPrecomputedTwoParticleForces;
    calculatePrecomputedTwoParticleForces();
}

void USCSIO2WaterPotential::initialize() {
    int highestAtomicNumber = +AtomTypes::Silicon;
    int numberOfAtomicNumbers = highestAtomicNumber+1;
    UnitConverter uc;

    m_activeAtomTypes.push_back(+AtomTypes::Oxygen);
    m_activeAtomTypes.push_back(+AtomTypes::Silicon);
    m_activeAtomTypes.push_back(+AtomTypes::Hydrogen);
    m_atomicNumberMap.resize(numberOfAtomicNumbers,-1);
    int atomMapIndex = 0; // First used atom will be the first index
    for(int atomType : m_activeAtomTypes) {
        m_atomicNumberMap.at(atomType) = atomMapIndex++;
    }


    m_configurationMap.resize(numberOfAtomicNumbers,vector<vector<int> >(numberOfAtomicNumbers, vector<int>(numberOfAtomicNumbers,0)));
    // Two particle configurations
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen).at(+AtomTypes::NoAtom) = +AtomConfiguration::O_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_O;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Silicon).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_Si;

    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Hydrogen).at(+AtomTypes::NoAtom) = +AtomConfiguration::O_H;
    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Oxygen).at(+AtomTypes::NoAtom) = +AtomConfiguration::O_H;

    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Hydrogen).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_H;
    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Silicon).at(+AtomTypes::NoAtom) = +AtomConfiguration::Si_H;

    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Hydrogen).at(+AtomTypes::NoAtom) = +AtomConfiguration::H_H;

    // Si-O-Si
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen) = +AtomConfiguration::Si_O_Si;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_Si;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_Si;


    // O-Si-O
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon) = +AtomConfiguration::O_Si_O;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen) = +AtomConfiguration::O_Si_O;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::Oxygen) = +AtomConfiguration::O_Si_O;

    // H_O_H
    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Hydrogen).at(+AtomTypes::Oxygen) = +AtomConfiguration::H_O_H;
    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Oxygen).at(+AtomTypes::Hydrogen) = +AtomConfiguration::H_O_H;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Hydrogen).at(+AtomTypes::Hydrogen) = +AtomConfiguration::H_O_H;

    // Si_O_H
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen).at(+AtomTypes::Hydrogen) = +AtomConfiguration::Si_O_H;
    m_configurationMap.at(+AtomTypes::Silicon).at(+AtomTypes::Hydrogen).at(+AtomTypes::Oxygen) = +AtomConfiguration::Si_O_H;

    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon).at(+AtomTypes::Hydrogen) = +AtomConfiguration::Si_O_H;
    m_configurationMap.at(+AtomTypes::Oxygen).at(+AtomTypes::Hydrogen).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_H;

    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Oxygen).at(+AtomTypes::Silicon) = +AtomConfiguration::Si_O_H;
    m_configurationMap.at(+AtomTypes::Hydrogen).at(+AtomTypes::Silicon).at(+AtomTypes::Oxygen) = +AtomConfiguration::Si_O_H;
}



void USCSIO2WaterPotential::calculatePrecomputedTwoParticleForces()
{
    m_waterParameters->createPrecomputedTables(m_numberOfPrecomputedTwoParticleForces, m_deltaR2);
    m_silicaParameters->createPrecomputedTables(m_numberOfPrecomputedTwoParticleForces, m_deltaR2);
    m_oneOverDeltaR2 = 1.0/m_deltaR2;
}

