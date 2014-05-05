#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <string>
#include <vector>

class AtomManager;
class Potential
{
private:
    std::string m_name;
protected:
    std::vector<double> m_systemLength;
    double m_potentialEnergy;
public:
    Potential();
    virtual void calculateForces(AtomManager &atomManager) = 0;
    std::string name() const;
    void setName(const std::string &name);
    std::vector<double> systemLength() const;
    void setSystemLength(const std::vector<double> &systemLength);
    double potentialEnergy() const;
    void setPotentialEnergy(double potentialEnergy);
};

#endif // POTENTIAL_H
