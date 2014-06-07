#include <statistics/totalenergysampler.h>
#include <statistics/kineticenergysampler.h>
#include <statistics/potentialenergysampler.h>
#include <includes.h>
using namespace Measurements;

StatisticalValue<double> TotalEnergySampler::value()
{
    return m_value;
}

shared_ptr<TotalEnergySampler> TotalEnergySampler::create(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<PotentialEnergySampler> potentialEnergySampler)
{
    return make_shared<TotalEnergySampler>(kineticEnergySampler, potentialEnergySampler);
}

TotalEnergySampler::TotalEnergySampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<PotentialEnergySampler> potentialEnergySampler) :
    m_kineticEnergySampler(kineticEnergySampler),
    m_potentialEnergySampler(potentialEnergySampler)
{
    addInput(kineticEnergySampler);
    addInput(potentialEnergySampler);
}

void TotalEnergySampler::action()
{
    double kineticEnergy = m_kineticEnergySampler.lock()->value().currentValueScalar();
    double potentialEnergy = m_potentialEnergySampler.lock()->value().currentValueScalar();
    double totalEnergy = kineticEnergy + potentialEnergy;

    m_value.addValue(totalEnergy);
}
