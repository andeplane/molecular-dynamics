#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
// private:
    public:
    void halfKick(System &system, double timestep);
    void move(System &system, double timestep);
    bool m_firstStep;

    void firstKick(System &system, const double &timestep);
    VelocityVerlet();
    virtual void integrate(System &system, const double &timestep);
};
