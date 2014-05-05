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

void VelocityVerlet::halfKick(System &system, double timestep)
{
    double dtHalf = timestep * 0.5;
    system.atomManager().atoms().iterate([&](Atom &atom) {
        atom.kick(dtHalf, atom.type()->massInverse());
    });
}

void VelocityVerlet::move(System &system, double timestep)
{
    system.atomManager().atoms().iterate([&](Atom &atom) {
        atom.move(timestep);
    });
}

void VelocityVerlet::firstKick(System &system, const double &timestep) {
    m_firstStep = false;

    system.atomManager().atoms().iterate([](Atom &atom) {
        atom.resetForce();
    });

    for(Potential *potential: system.potentials()) {
        potential->calculateForces(system.atomManager());
    }

    halfKick(system,timestep);
}

void VelocityVerlet::integrate(System &system, const double &timestep)
{
    if(m_firstStep) firstKick(system, timestep);

    halfKick(system, timestep);
    move(system, timestep);
    system.topology().MPIMove(system);

    system.atomManager().atoms().iterate([](Atom &atom) {
        atom.resetForce();
    });

    for(Potential *potential: system.potentials()) {
        potential->calculateForces(system.atomManager());
    }

    halfKick(system, timestep);
}
