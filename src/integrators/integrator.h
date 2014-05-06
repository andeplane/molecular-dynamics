#pragma once
enum class Integrators {VelocityVerlet = 0};

class System;
class Integrator
{
public:
    Integrator();
    virtual void integrate(System &system, const double &timestep) = 0;
};
