#include <potentials/lennardjonespotential.h>
#include <system.h>
#include <functional>
#include <vector>
#include <unitconverter.h>

unsigned long LennardJonesPotential::numberOfComputedPairs() const
{
    return m_numberOfComputedPairs;
}


unsigned long LennardJonesPotential::numberOfComputedPairsWithinCutoffDistance() const
{
    return m_numberOfComputedPairsWithinCutoffDistance;
}

bool LennardJonesPotential::calculateForcesBetweenAllPairsWithMinimumImageConvention() const
{
    return m_calculateForcesBetweenAllPairsWithMinimumImageConvention;
}

void LennardJonesPotential::setCalculateForcesBetweenAllPairsWithMinimumImageConvention(bool calculateForcesBetweenAllPairsWithMinimumImageConvention)
{
    m_calculateForcesBetweenAllPairsWithMinimumImageConvention = calculateForcesBetweenAllPairsWithMinimumImageConvention;
}

LennardJonesPotential::LennardJonesPotential() :
    m_sigma(UnitConverter::lengthFromAngstroms(3.405)),
    m_epsilon(UnitConverter::energyFromSI(1.65088e-21)),
    m_numberOfComputedPairs(0),
    m_numberOfComputedPairsWithinCutoffDistance(0),
    m_calculateForcesBetweenAllPairsWithMinimumImageConvention(false),
    m_potentialEnergyCorrection(0)
{
    setName("Lennard Jones");
    using namespace std::placeholders;

    m_iteratorAllPairs.setTwoParticleAction(std::bind(&LennardJonesPotential::twoParticleActionMinimumImageConvention, this, _1, _2));
    m_iteratorDefault.setTwoParticleAction(std::bind(&LennardJonesPotential::twoParticleAction, this, _1, _2));
    setCutoffDistance(2.5*UnitConverter::lengthFromAngstroms(3.405));
}

void LennardJonesPotential::twoParticleActionMinimumImageConvention(Atom *atom1, Atom *atom2) {
    if(m_systemLength.x() == 0) {
        std::cout << "Warning, calculating forces with Lennard Jones potential using minimum image convention systemLength set to zero, aborting." << std::endl;
        return;
    }
    m_numberOfComputedPairs++;

    double dr[3];
    double cutoffDistanceSquared = m_cutoffDistance*m_cutoffDistance;
    dr[0] = atom1->position[0] - atom2->position[0];
    dr[1] = atom1->position[1] - atom2->position[1];
    dr[2] = atom1->position[2] - atom2->position[2];

    if(dr[0] > m_systemLength[0]*0.5) dr[0] -= m_systemLength[0];
    if(dr[0] < -m_systemLength[0]*0.5) dr[0] += m_systemLength[0];

    if(dr[1] > m_systemLength[1]*0.5) dr[1] -= m_systemLength[1];
    if(dr[1] < -m_systemLength[1]*0.5) dr[1] += m_systemLength[1];

    if(dr[2] > m_systemLength[2]*0.5) dr[2] -= m_systemLength[2];
    if(dr[2] < -m_systemLength[2]*0.5) dr[2] += m_systemLength[2];

    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
#ifdef DEBUG
    if(dr2 == 0) {
        std::cout << "Error in LennardJonesPotential::twoParticleAction. Relative distance is zero:" << std::endl;
        std::cout << *atom1 << endl;
        std::cout << *atom2 << endl;
        return;
    }
#endif

    if (dr2<cutoffDistanceSquared) {
        m_numberOfComputedPairsWithinCutoffDistance++;
        double dr2Inverse = 1.0/dr2;
        double sigmaOverDr2 = m_sigma*m_sigma*dr2Inverse;
        double sigmaOverDr6 = sigmaOverDr2*sigmaOverDr2*sigmaOverDr2;
        double force = (2*sigmaOverDr6-1)*sigmaOverDr6*dr2Inverse*24*m_epsilon;

        atom1->force[0] += force*dr[0];
        atom1->force[1] += force*dr[1];
        atom1->force[2] += force*dr[2];
        atom2->force[0] -= force*dr[0];
        atom2->force[1] -= force*dr[1];
        atom2->force[2] -= force*dr[2];
    }
}

void LennardJonesPotential::twoParticleAction(Atom *atom1, Atom *atom2)
{
    m_numberOfComputedPairs++;

    double dr[3];
    double cutoffDistanceSquared = m_cutoffDistance*m_cutoffDistance;
    dr[0] = atom1->position[0] - atom2->position[0];
    dr[1] = atom1->position[1] - atom2->position[1];
    dr[2] = atom1->position[2] - atom2->position[2];
    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    if(dr2 < 1e-5) {
        std::cout << "Error in LennardJonesPotential::twoParticleAction. Relative distance is zero:" << std::endl;
        std::cout << *atom1 << endl;
        std::cout << *atom2 << endl;
        return;
    }

    if (dr2<cutoffDistanceSquared) {
        m_numberOfComputedPairsWithinCutoffDistance++;
        double dr2Inverse = 1.0/dr2;
        double sigmaOverDr2 = m_sigma*m_sigma*dr2Inverse;
        double sigmaOverDr6 = sigmaOverDr2*sigmaOverDr2*sigmaOverDr2;
        double force = (2*sigmaOverDr6-1)*sigmaOverDr6*dr2Inverse*24*m_epsilon;

        int numberOfGhosts = atom1->ghost() + atom2->ghost();
        if(numberOfGhosts == 0) {
            m_potentialEnergy += 4*m_epsilon*sigmaOverDr6*(sigmaOverDr6 - 1) - m_potentialEnergyCorrection;
        } else {
            m_potentialEnergy += 0.5*(4*m_epsilon*sigmaOverDr6*(sigmaOverDr6 - 1) - m_potentialEnergyCorrection);
        }

        atom1->force[0] += force*dr[0];
        atom1->force[1] += force*dr[1];
        atom1->force[2] += force*dr[2];
        atom2->force[0] -= force*dr[0];
        atom2->force[1] -= force*dr[1];
        atom2->force[2] -= force*dr[2];
    }
}

double LennardJonesPotential::cutoffDistance() const
{
    return m_cutoffDistance;
}

void LennardJonesPotential::setCutoffDistance(double cutoffDistance)
{
    m_cutoffDistance = cutoffDistance;

    // Remember to shift the potential energy so it is zero at the cutoff distance
    m_potentialEnergyCorrection = 4*m_epsilon*pow(m_sigma/m_cutoffDistance,6.0)*(pow(m_sigma/m_cutoffDistance,6.0) - 1);
}

double LennardJonesPotential::sigma() const
{
    return m_sigma;
}

void LennardJonesPotential::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJonesPotential::epsilon() const
{
    return m_epsilon;
}

void LennardJonesPotential::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJonesPotential::calculateForces(AtomManager &atomManager)
{
    setPotentialEnergy(0.0);
    atomManager.setCutoffDistance(m_cutoffDistance);

    if(m_calculateForcesBetweenAllPairsWithMinimumImageConvention) {
        m_iteratorAllPairs.setLoopThroughGhosts(false);
        m_iteratorAllPairs.iterate(atomManager);
    } else {
        m_iteratorDefault.iterate(atomManager);
    }
}
