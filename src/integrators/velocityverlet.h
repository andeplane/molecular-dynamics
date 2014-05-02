#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System &system, double timestep);
    void move(System &system, double timestep);
    bool m_firstStep;

    void firstKick(System &system, const double &timestep);
public:
    VelocityVerlet();
    virtual void integrate(System &system, const double &timestep);
};

#endif // VELOCITYVERLET_H
