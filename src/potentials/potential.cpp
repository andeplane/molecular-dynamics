#include <potentials/potential.h>


std::string Potential::name() const
{
    return m_name;
}

void Potential::setName(const std::string &name)
{
    m_name = name;
}
Potential::Potential()
{
    m_name = "Unnamed potential";
}
