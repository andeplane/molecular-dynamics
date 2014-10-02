#pragma once
#include <modifiers/modifier.h>
#include <statistics/temperaturesampler.h>
namespace Modifiers {
class BerendsenThermostat : public Modifier
{
private:
    double m_deltaT;
    double m_relaxationTime;
    double m_temperature;
    StatisticalValue<double> m_value;
    weak_ptr<Measurements::TemperatureSampler> m_temperatureSampler;
    weak_ptr<System> m_system;
public:
    static shared_ptr<BerendsenThermostat> create(double temperature, double deltaT, double relaxationTime, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system);
    BerendsenThermostat(double temperature, double deltaT, double relaxationTime, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system);
    virtual void action();
    StatisticalValue<double> value();
};
}
