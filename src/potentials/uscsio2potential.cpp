#include "uscsio2potential.h"
#include <atomtype.h>

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

vector<double> USCSIO2Potential::cutoffDistances() const
{
    return m_cutoffDistances;
}


double USCSIO2Potential::maxCutoffDistance() const
{
    return m_maxCutoffDistance;
}

void USCSIO2Potential::initialize() {
    m_atomicNumberMap.resize(30,0);
    m_atomicNumberMap.at(int(AtomTypes::Oxygen)) = 0;
    m_atomicNumberMap.at(int(AtomTypes::Silicon)) = 1;
}

USCSIO2Potential::USCSIO2Potential()
{
    setName("USC SiO2");
    using namespace std::placeholders;

    m_iteratorDefault.setTwoParticleAction(std::bind(&USCSIO2Potential::twoParticleAction, this, _1, _2));
    m_iteratorDefault.setThreeParticleAction(std::bind(&USCSIO2Potential::threeParticleAction, this, _1, _2, _3));

}
