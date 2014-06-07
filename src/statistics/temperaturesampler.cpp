#include <statistics/temperaturesampler.h>
#include <statistics/kineticenergysampler.h>


StatisticalValue<double> TemperatureSampler::value()
{
    return m_value;
}

TemperatureSampler::TemperatureSampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system) :
    m_kineticEnergySampler(kineticEnergySampler),
    m_system(system)
{
    addDependency(kineticEnergySampler);
}

void TemperatureSampler::action()
{
    double kineticEnergy = m_kineticEnergySampler.lock()->value().currentValueScalar();
    double temperature = 2.0/(3.0*m_system.lock()->numberOfAtoms())*kineticEnergy;
    m_value.addValue(temperature);
}
