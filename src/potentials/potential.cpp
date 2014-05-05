#include <potentials/potential.h>


std::string Potential::name() const
{
    return m_name;
}

void Potential::setName(const std::string &name)
{
    m_name = name;
}

std::vector<double> Potential::systemLength() const
{
    return m_systemLength;
}

void Potential::setSystemLength(const std::vector<double> &systemLength)
{
    m_systemLength = systemLength;
}

double Potential::potentialEnergy() const
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}
Potential::Potential():
    m_name("Unnamed potential"),
    m_potentialEnergy(0)
{
    m_systemLength.resize(3,0);
}
