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
// #define DEBUGVELOCITYVERLET

#ifdef DEBUGVELOCITYVERLET
void VelocityVerlet::integrate(System &system, const double &timestep)
{
    cout << "Beginning of timestep" << endl;
    cout << system.atomManager() << endl;
    halfKick(system, timestep);
    cout << "After half kick" << endl;
    cout << system.atomManager() << endl;

    move(system, timestep);
    cout << "After move" << endl;
    cout << system.atomManager() << endl;

    system.topology().MPIMove(system);
    cout << "After MPI Move" << endl;
    cout << system.atomManager() << endl;

    system.atomManager().atoms().iterate([](Atom &atom) {
        atom.resetForce();
    });

    cout << "After reset force" << endl;
    cout << system.atomManager() << endl;

    for(Potential *potential: system.potentials()) {
        potential->calculateForces(system.atomManager());
    }

    cout << "After calculating forces" << endl;
    cout << system.atomManager() << endl;

    halfKick(system, timestep);
    cout << "After half kick" << endl;
    cout << system.atomManager() << endl;
}
#else

void VelocityVerlet::firstKick(System &system, const double &timestep) {
    m_firstStep = false;
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
#endif
