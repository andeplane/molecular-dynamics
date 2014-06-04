#pragma once
#include <statistics/statisticalproperty.h>

class KineticEnergySampler : public StatisticalProperty
{
private:
    StatisticalValue<double> m_value;
    shared_ptr<System> m_system;
public:
    KineticEnergySampler(shared_ptr<System> system);
    StatisticalValue<double> value();
    virtual void action();
};
