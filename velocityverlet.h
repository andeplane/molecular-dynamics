#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include <integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System &system, double timestep);
    void move(System &system, double timestep);

public:
    VelocityVerlet();

};

#endif // VELOCITYVERLET_H
