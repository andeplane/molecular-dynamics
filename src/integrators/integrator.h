#ifndef INTEGRATOR_H
#define INTEGRATOR_H

enum class Integrators {VelocityVerlet = 0};

class System;
class Integrator
{
public:
    Integrator();
    virtual void integrate(System &system, const double &timestep) = 0;
};

#endif // INTEGRATOR_H
