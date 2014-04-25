#include "lennardjonespotential.h"

LennardJonesPotential::LennardJonesPotential() :
    m_sigma(1),
    m_cutoffRadius(3.0),
    m_epsilon(1)
{

}

double LennardJonesPotential::cutoffRadius() const
{
    return m_cutoffRadius;
}

void LennardJonesPotential::setCutoffRadius(double cutoffRadius)
{
    m_cutoffRadius = cutoffRadius;
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

void LennardJonesPotential::calculateForces()
{

}
