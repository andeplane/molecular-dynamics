#include <potentials/lennardjonespotential.h>
#include <system.h>


double LennardJonesPotential::cutoffDistance() const
{
    return m_cutoffDistance;
}

void LennardJonesPotential::setCutoffDistance(double cutoffDistance)
{
    m_cutoffDistance = cutoffDistance;
}
LennardJonesPotential::LennardJonesPotential() :
    m_sigma(1),
    m_epsilon(1),
    m_cutoffDistance(2.5)
{

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
#include <iostream>
using namespace std;
void LennardJonesPotential::calculateForces(AtomManager &atomManager)
{
    double cutoffDistanceSquared = m_cutoffDistance*m_cutoffDistance;

    m_iterator.setTwoParticleAction([&](Atom *atom1, Atom *atom2) {
        double dr[3];
        dr[0] = atom1->position[0] - atom2->position[0];
        dr[1] = atom1->position[1] - atom2->position[1];
        dr[2] = atom1->position[2] - atom2->position[2];
        double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

        if (dr2<cutoffDistanceSquared) {
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
    });

    m_iterator.iterate(atomManager);
}
