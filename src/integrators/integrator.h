#pragma once
#include <memory>
using std::shared_ptr;

enum class Integrators {VelocityVerlet = 0};

class System;
class Integrator
{
public:
    Integrator();
    virtual void integrate(shared_ptr<System> system, const double &timestep) = 0;
};
