#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
// private:
    public:
    void halfKick(shared_ptr<System> system, double timestep);
    void move(shared_ptr<System> system, double timestep);
    bool m_firstStep;

    void firstKick(shared_ptr<System> system, const double &timestep);
    VelocityVerlet();
    virtual void integrate(shared_ptr<System> system, const double &timestep);
};
