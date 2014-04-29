#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System &system, double timestep);
    void move(System &system, double timestep);

public:
    VelocityVerlet();
    virtual void integrate(System &system, double timestep);
};

#endif // VELOCITYVERLET_H
