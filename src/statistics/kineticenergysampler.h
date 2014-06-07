#pragma once
#include <statistics/statisticalproperty.h>

class KineticEnergySampler : public StatisticalProperty
{
private:
    StatisticalValue<double> m_value;
    weak_ptr<System> m_system;
public:
    KineticEnergySampler(weak_ptr<System> system);
    StatisticalValue<double> value();
    virtual void action();
};
