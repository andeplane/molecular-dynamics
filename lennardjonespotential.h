#ifndef LENNARDJONESPOTENTIAL_H
#define LENNARDJONESPOTENTIAL_H
#include <potential.h>

class LennardJonesPotential : public Potential
{
private:
    double m_cutoffRadius;
    double m_sigma;
    double m_epsilon;
public:
    LennardJonesPotential();
    double cutoffRadius() const;
    void setCutoffRadius(double cutoffRadius);
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
    virtual void calculateForces();
};

#endif // LENNARDJONESPOTENTIAL_H
