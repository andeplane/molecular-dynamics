#pragma once
#include <string>
#include <vector>

enum class AtomConfiguration {NotUsed = 0, Si_O=1, Si_Si=2, O_O=3, O_H=4, H_H=5, Si_H=6, O_Si_O=7, Si_O_Si=8, H_O_H=9, Si_O_H=10, NumberOfConfigurations=11};

inline int operator + (AtomConfiguration val) {
    return static_cast<int>(val);
}

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
