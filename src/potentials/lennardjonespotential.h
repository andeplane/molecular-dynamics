#ifndef LENNARDJONESPOTENTIAL_H
#define LENNARDJONESPOTENTIAL_H
#include <potentials/potential.h>

class LennardJonesPotential : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
public:
    LennardJonesPotential();
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
    virtual void calculateForces(System &system);
};

#endif // LENNARDJONESPOTENTIAL_H
