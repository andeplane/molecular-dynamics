#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class System;

class Integrator
{
public:
    Integrator();
    virtual void integrate(System &system, double timestep) = 0;
};

#endif // INTEGRATOR_H
