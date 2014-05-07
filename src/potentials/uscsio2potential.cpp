/* USC SiO2 potential
 * [1] Vashishta, P., et al. "Interaction potential for SiO 2: a molecular-dynamics study of structural correlations." Physical Review B 41.17 (1990): 12197.
 */

#include <potentials/uscsio2potential.h>
#include <unitconverter.h>

int operator + (AtomConfiguration val) {
    return static_cast<int>(val);
}

void USCSIO2Potential::calculateForces(AtomManager &atomManager)
{
    atomManager.setCutoffDistance(m_maxCutoffDistance);
    m_iteratorDefault.iterate(atomManager);
}

void USCSIO2Potential::twoParticleAction(Atom *atom1, Atom *atom2)
{
    double dr[3];

    dr[0] = atom1->position[0] - atom2->position[0];
    dr[1] = atom1->position[1] - atom2->position[1];
    dr[2] = atom1->position[2] - atom2->position[2];

    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

#ifdef DEBUG
    if(dr2 < 1e-5) {
        std::cout << "Error in USCSIO2Potential::twoParticleAction. Relative distance is near zero:" << std::endl;
        std::cout << *atom1 << endl;
        std::cout << *atom2 << endl;
        return;
    }
#endif
}

void USCSIO2Potential::threeParticleAction(Atom *atom1, Atom *atom2, Atom *atom3)
{

}


double USCSIO2Potential::maxCutoffDistance() const
{
    return m_maxCutoffDistance;
}

void USCSIO2Potential::initialize() {
    UnitConverter uc;

    vector<int> activeAtomTypes;
    activeAtomTypes.push_back(+AtomTypes::Oxygen);
    activeAtomTypes.push_back(+AtomTypes::Silicon);
    m_atomicNumberMap.resize(30,-1);
    int atomMapIndex = 0; // First used atom will be the first index
    for(int atomType : activeAtomTypes) {
        m_atomicNumberMap.at(atomType) = atomMapIndex++;
    }

    B_ijk.resize(+AtomConfiguration::NumberOfConfigurations,0);
    theta_zero.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r_0.resize(+AtomConfiguration::NumberOfConfigurations,0);
    ksi.resize(+AtomConfiguration::NumberOfConfigurations,0);
    C_ijk.resize(+AtomConfiguration::NumberOfConfigurations,0);
    eta_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    H_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    D_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    W_ij.resize(+AtomConfiguration::NumberOfConfigurations,0);
    cutoffDistances.resize(+AtomConfiguration::NumberOfConfigurations,0);
    Z_i.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r1s.resize(+AtomConfiguration::NumberOfConfigurations,0);
    r4s.resize(+AtomConfiguration::NumberOfConfigurations,0);

    // ALL THESE VALUES ARE FROM TABLE 1 IN [1]
    r1s.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(4.43);
    r4s.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(2.5);

    Z_i.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(2.5);

    eta_ij.at(+AtomConfiguration::Si_O) = 9;
    eta_ij.at(+AtomConfiguration::Si_Si) = 11;
    eta_ij.at(+AtomConfiguration::O_O) = 7;

    cutoffDistances.at(+AtomConfiguration::Si_O) = uc.lengthFromAngstroms(5.5);
    cutoffDistances.at(+AtomConfiguration::Si_Si) = uc.lengthFromAngstroms(5.5);
    cutoffDistances.at(+AtomConfiguration::O_O) = uc.lengthFromAngstroms(5.5);

    D_ij.at(+AtomConfiguration::Si_O) = 3.456*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),4.0);
    D_ij.at(+AtomConfiguration::Si_Si) = 0.0;
    D_ij.at(+AtomConfiguration::O_O) = 1.728*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),4.0);

    H_ij.at(+AtomConfiguration::Si_O) = 78.3143*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_O));
    H_ij.at(+AtomConfiguration::Si_Si) = 0.39246*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::Si_Si));
    H_ij.at(+AtomConfiguration::O_O) = 355.5263*uc.energyFromEv(1.0)*pow(uc.lengthFromAngstroms(1.0),eta_ij.at(+AtomConfiguration::O_O));

    B_ijk.at(+AtomConfiguration::O_Si_O) = uc.energyFromEv(4.993);
    B_ijk.at(+AtomConfiguration::Si_O_Si) = uc.energyFromEv(19.972);

    theta_zero.at(+AtomConfiguration::O_Si_O) = uc.degreesToRadians(109.47);
    theta_zero.at(+AtomConfiguration::Si_O_Si) = uc.degreesToRadians(141.0);

    r_0.at(+AtomConfiguration::O_Si_O) = uc.lengthFromAngstroms(2.60);
    r_0.at(+AtomConfiguration::Si_O_Si) = uc.lengthFromAngstroms(2.60);

    ksi.at(+AtomConfiguration::O_Si_O) = uc.lengthFromAngstroms(1.0);
    ksi.at(+AtomConfiguration::Si_O_Si) = uc.lengthFromAngstroms(1.0);

    C_ijk.at(+AtomConfiguration::O_Si_O) = 0;
    C_ijk.at(+AtomConfiguration::Si_O_Si) = 0;
}

USCSIO2Potential::USCSIO2Potential()
{
    setName("USC SiO2");
    using namespace std::placeholders;

    m_iteratorDefault.setTwoParticleAction(std::bind(&USCSIO2Potential::twoParticleAction, this, _1, _2));
    m_iteratorDefault.setThreeParticleAction(std::bind(&USCSIO2Potential::threeParticleAction, this, _1, _2, _3));

}
