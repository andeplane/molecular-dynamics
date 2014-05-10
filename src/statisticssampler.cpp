#include <statisticssampler.h>
#include <system.h>
#include <atom.h>
#include <atomtype.h>
#include <potentials/potential.h>

StatisticsSampler::StatisticsSampler()
{

}

double StatisticsSampler::calculateKineticEnergy(System &system)
{
    double kineticEnergy = 0;
    system.atomManager().atoms().iterate([&](Atom &atom) {
        kineticEnergy += 0.5*atom.type()->mass()*(atom.velocity[0]*atom.velocity[0] + atom.velocity[1]*atom.velocity[1] + atom.velocity[2]*atom.velocity[2]);
    });

    return kineticEnergy;
}

double StatisticsSampler::calculatePotentialEnergy(System &system)
{
    double potentialEnergy = 0;
    for(Potential *potential : system.potentials()) {
        potentialEnergy += potential->potentialEnergy();
    }

    return potentialEnergy;
}

double StatisticsSampler::calculateTotalEnergy(System &system)
{
    return calculatePotentialEnergy(system) + calculateKineticEnergy(system);
}

vector<double> StatisticsSampler::calculateTotalMomentum(System &system)
{
    vector<double> momentum(3,0);
    system.atomManager().atoms().iterate([&](Atom &atom) {
        momentum[0] += atom.velocity[0]*atom.type()->mass();
        momentum[1] += atom.velocity[1]*atom.type()->mass();
        momentum[2] += atom.velocity[2]*atom.type()->mass();
    });

    return momentum;
}
