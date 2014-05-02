#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <string>

class AtomManager;
class Potential
{
private:
    std::string m_name;
public:
    Potential();
    virtual void calculateForces(AtomManager &atomManager) = 0;
    std::string name() const;
    void setName(const std::string &name);
};

#endif // POTENTIAL_H
