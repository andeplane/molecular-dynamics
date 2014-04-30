#include <integrators/velocityverlet.h>
#include <system.h>
#include <topology.h>
#include <potentials/potential.h>
#include <atom.h>
#include <iostream>
#include <atommanager.h>

using namespace std;

VelocityVerlet::VelocityVerlet()
{

}

void VelocityVerlet::halfKick(System &system, double timestep)
{
//    double dtHalf = timestep * 0.5;
//    for(Atom *atom : system.atoms()) {
//        atom->kick(dtHalf, atom->type()->massInverse());
//    }
}

void VelocityVerlet::move(System &system, double timestep)
{

}

void VelocityVerlet::integrate(System &system, const double &timestep)
{
    halfKick(system, timestep);
    move(system, timestep);
    system.topology().MPIMove(system);
    system.topology().MPICopy(system,system.atomManager().cutoffDistance());
    system.resetForces();

    for(Potential *potential: system.potentials()) {
        potential->calculateForces(system.atomManager());
    }

    halfKick(system, timestep);
}
