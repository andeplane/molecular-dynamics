#include "velocityverlet.h"
#include <system.h>
#include <topology.h>
#include <potential.h>
#include <atom.h>

VelocityVerlet::VelocityVerlet()
{

}

void VelocityVerlet::halfKick(System &system, double timestep)
{
    double dtHalf = timestep * 0.5;
    double oneOverMass = 1.0/39.948;
    for(Atom atom : system.atoms()) {
        atom.kick(dtHalf, oneOverMass);
    }
}

void VelocityVerlet::move(System &system, double timestep)
{
    for(Atom atom : system.atoms()) {
        atom.move(timestep);
    }
}

void VelocityVerlet::integrate(System &system, double timestep) {
    halfKick(system, timestep);
    move(system, timestep);
    system.topology().MPIMove(system);
    system.topology().MPICopy(system);
    system.resetForces();

    for(Potential potential: system.potentials()) {
        potential.calculateForces(system.atoms());
    }

    halfKick(system, timestep);
}
