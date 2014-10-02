#include <modifiers/berendsenthermostat.h>
#include <particles/atom.h>

using namespace Modifiers;

void BerendsenThermostat::action() {
    double temperature = m_temperatureSampler.lock().get()->value().currentValueScalar();
    double factor = sqrt(1 + m_relaxationTime/m_deltaT*(m_temperature/temperature - 1));
    m_system.lock().get()->atomManager().atoms().iterate([&](Atom &atom) {
        atom.velocity[0] *= factor;
        atom.velocity[1] *= factor;
        atom.velocity[2] *= factor;
    });
}



shared_ptr<BerendsenThermostat> BerendsenThermostat::create(double temperature, double deltaT, double relaxationTime, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system)
{
    return std::make_shared<BerendsenThermostat>(temperature, deltaT, relaxationTime, temperatureSampler, system);
}

BerendsenThermostat::BerendsenThermostat(double temperature, double deltaT, double relaxationTime, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system) :
    m_deltaT(deltaT),
    m_relaxationTime(relaxationTime),
    m_temperature(temperature),
    m_temperatureSampler(temperatureSampler),
    m_system(system)
{
    addInput(temperatureSampler);
}

StatisticalValue<double> BerendsenThermostat::value()
{
    return m_value;
}
