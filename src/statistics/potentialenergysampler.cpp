#include <statistics/potentialenergysampler.h>
#include <potentials/potential.h>


StatisticalValue<double> PotentialEnergySampler::value()
{
    return m_value;
}

PotentialEnergySampler::PotentialEnergySampler(shared_ptr<System> system) :
    m_system(system)
{

}

void PotentialEnergySampler::action()
{
    double potentialEnergy = 0;
    for(Potential *potential : m_system->potentials()) {
        potentialEnergy += potential->potentialEnergy();
    }

    m_value.addValue(potentialEnergy);
}
