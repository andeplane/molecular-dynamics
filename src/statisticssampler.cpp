#include <statisticssampler.h>
#include <system.h>
#include <particles/atom.h>
#include <atomtype.h>
#include <potentials/potential.h>

StatisticsSampler::StatisticsSampler()
{

}

double StatisticsSampler::calculateKineticEnergy(shared_ptr<System> system)
{
    double kineticEnergy = 0;
    system.get()->atomManager().atoms().iterate([&](Atom &atom) {
        kineticEnergy += 0.5*atom.type()->mass()*(atom.velocity[0]*atom.velocity[0] + atom.velocity[1]*atom.velocity[1] + atom.velocity[2]*atom.velocity[2]);
    });

    return kineticEnergy;
}

double StatisticsSampler::calculatePotentialEnergy(shared_ptr<System> system)
{
    double potentialEnergy = 0;

    for(Potential *potential : system.get()->potentials()) {
        potentialEnergy += potential->potentialEnergy();
    }

    return potentialEnergy;
}

double StatisticsSampler::calculateTemperature(shared_ptr<System> system)
{
    double kineticEnergy = calculateKineticEnergy(system);

    double temperature = 2.0/(3*system.get()->numberOfAtoms())*kineticEnergy;

    return temperature;
}

double StatisticsSampler::calculateTotalEnergy(shared_ptr<System> system)
{
    return calculatePotentialEnergy(system) + calculateKineticEnergy(system);
}

vector<double> StatisticsSampler::calculateTotalMomentum(shared_ptr<System> system)
{
    vector<double> momentum(3,0);
    system.get()->atomManager().atoms().iterate([&](Atom &atom) {
        momentum[0] += atom.velocity[0]*atom.type()->mass();
        momentum[1] += atom.velocity[1]*atom.type()->mass();
        momentum[2] += atom.velocity[2]*atom.type()->mass();
    });

    return momentum;
}
