#include <statistics/potentialenergysampler.h>
#include <potentials/potential.h>
#include <includes.h>

StatisticalValue<double> Measurements::PotentialEnergySampler::value()
{
    return m_value;
}

shared_ptr<Measurements::PotentialEnergySampler> Measurements::PotentialEnergySampler::create(shared_ptr<System> system)
{
    return make_shared<Measurements::PotentialEnergySampler>(system);
}

Measurements::PotentialEnergySampler::PotentialEnergySampler(shared_ptr<System> system) :
    m_system(system)
{

}

void Measurements::PotentialEnergySampler::action()
{
    double potentialEnergy = 0;
    for(Potential *potential : m_system.lock()->potentials()) {
        potentialEnergy += potential->potentialEnergy();
    }

    m_value.addValue(potentialEnergy);
}
