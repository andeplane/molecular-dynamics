#include <statistics/temperaturesampler.h>
#include <statistics/kineticenergysampler.h>
#include <includes.h>
using namespace Measurements;

StatisticalValue<double> TemperatureSampler::value()
{
    return m_value;
}

shared_ptr<TemperatureSampler> TemperatureSampler::create(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system)
{
    return make_shared<TemperatureSampler>(kineticEnergySampler, system);
}

TemperatureSampler::TemperatureSampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system) :
    m_kineticEnergySampler(kineticEnergySampler),
    m_system(system)
{
    addInput(kineticEnergySampler);
}

void TemperatureSampler::action()
{
    double kineticEnergy = m_kineticEnergySampler.lock()->value().currentValueScalar();
    double temperature = 2.0/(3.0*m_system.lock()->numberOfAtoms())*kineticEnergy;
    m_value.addValue(temperature);
}
