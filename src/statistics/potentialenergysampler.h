#pragma once
#include <statistics/statisticalproperty.h>
class System;

class PotentialEnergySampler : public StatisticalProperty
{
private:
    weak_ptr<System> m_system;
    StatisticalValue<double> m_value;

public:
    PotentialEnergySampler(weak_ptr<System> system);
    virtual void action();
    StatisticalValue<double> value();
};
