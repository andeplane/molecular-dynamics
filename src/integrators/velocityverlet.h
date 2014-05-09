#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System &system, double timestep);
    void move(System &system, double timestep);
    bool m_firstStep;
public:

    void firstKick(System &system, const double &timestep);
    VelocityVerlet();
    virtual void integrate(System &system, const double &timestep);
};
