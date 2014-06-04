#pragma once
#include <memory>
#include <statistics/statisticalproperty.h>
using std::shared_ptr;

class System;
class KineticEnergySampler;

class TemperatureSampler : public StatisticalProperty
{
private:
    StatisticalValue<double> m_value;
    shared_ptr<KineticEnergySampler> m_kineticEnergySampler;
    shared_ptr<System> m_system;
public:
    TemperatureSampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system);
    virtual void action();
    StatisticalValue<double> value();
};
