#include <potentials/potential.h>
using CompPhys::vec3;

std::string Potential::name() const
{
    return m_name;
}

void Potential::setName(const std::string &name)
{
    m_name = name;
}

vec3 Potential::systemLength() const
{
    return m_systemLength;
}

void Potential::setSystemLength(CompPhys::vec3 systemLength)
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

}
