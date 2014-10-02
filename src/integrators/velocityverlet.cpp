#include <integrators/velocityverlet.h>
#include <system.h>
#include <topology.h>
#include <potentials/potential.h>
#include <atom.h>
#include <iostream>
#include <atommanager.h>

using namespace std;

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true)
{

}

void VelocityVerlet::halfKick(shared_ptr<System> system, double timestep)
{
    double dtHalf = timestep * 0.5;
    system->atomManager().atoms().iterate([&](Atom &atom) {
        atom.kick(dtHalf, atom.type()->massInverse());
    });
}

void VelocityVerlet::move(shared_ptr<System> system, double timestep)
{
    system->atomManager().atoms().iterate([&](Atom &atom) {
        atom.move(timestep);
    });
}

void VelocityVerlet::firstKick(shared_ptr<System> system, const double &timestep) {
    m_firstStep = false;

    system->atomManager().atoms().iterate([](Atom &atom) {
        atom.resetForce();
    });

    for(Potential *potential: system->potentials()) {
        potential->calculateForces(system->atomManager());
    }

    halfKick(system,timestep);
}

void VelocityVerlet::integrate(shared_ptr<System> system, const double &timestep)
{
    if(m_firstStep) firstKick(system, timestep);
    else halfKick(system, timestep);
    move(system, timestep);
    system->topology().MPIMove(system->atomManager());

    system->atomManager().atoms().iterate([](Atom &atom) {
        atom.resetForce();
    });

    for(Potential *potential: system->potentials()) {
        potential->calculateForces(system->atomManager());
    }

    halfKick(system, timestep);
}
