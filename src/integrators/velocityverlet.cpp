#include <integrators/velocityverlet.h>
#include <system.h>
#include <topology.h>
#include <potentials/potential.h>
#include <atom.h>
#include <iostream>
using namespace std;

VelocityVerlet::VelocityVerlet()
{

}

void VelocityVerlet::halfKick(System &system, double timestep)
{
    double dtHalf = timestep * 0.5;
    double oneOverMass = 1.0/39.948;
    for(Atom *atom : system.atoms()) {
        atom->kick(dtHalf, atom->type()->massInverse());
    }
}

void VelocityVerlet::move(System &system, double timestep)
{
    for(Atom *atom : system.atoms()) {
        atom->move(timestep);
    }
}

void VelocityVerlet::integrate(System &system, double timestep)
{
    halfKick(system, timestep);
    move(system, timestep);
    system.topology().MPIMove(system);
    system.topology().MPICopy(system,system.cutoffDistance());
    system.resetForces();

    for(Potential *potential: system.potentials()) {
        potential->calculateForces(system);
    }

    halfKick(system, timestep);
}
