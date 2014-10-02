#pragma once
#include <string>
#include <vector>
#include <utils/vec3.h>

class AtomManager;
class Potential
{
private:
    std::string m_name;
protected:
    CompPhys::vec3 m_systemLength;
    double m_potentialEnergy;
public:
    Potential();
    virtual void calculateForces(AtomManager &atomManager) = 0;
    std::string name() const;
    void setName(const std::string &name);
    CompPhys::vec3 systemLength() const;
    void setSystemLength(CompPhys::vec3 systemLength);
    double potentialEnergy() const;
    void setPotentialEnergy(double potentialEnergy);
};
