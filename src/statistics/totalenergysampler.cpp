#include <statistics/totalenergysampler.h>
#include <statistics/kineticenergysampler.h>
#include <statistics/potentialenergysampler.h>

StatisticalValue<double> TotalEnergySampler::value()
{
    return m_value;
}

TotalEnergySampler::TotalEnergySampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<PotentialEnergySampler> potentialEnergySampler) :
    m_kineticEnergySampler(kineticEnergySampler),
    m_potentialEnergySampler(potentialEnergySampler)
{
    addDependency(kineticEnergySampler);
    addDependency(potentialEnergySampler);
}

void TotalEnergySampler::action()
{
    double kineticEnergy = m_kineticEnergySampler.lock()->value().currentValueScalar();
    double potentialEnergy = m_potentialEnergySampler.lock()->value().currentValueScalar();
    double totalEnergy = kineticEnergy + potentialEnergy;

    m_value.addValue(totalEnergy);
}
