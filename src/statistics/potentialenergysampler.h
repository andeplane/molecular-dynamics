#pragma once
#include <statistics/statisticalproperty.h>
class System;

class PotentialEnergySampler : public StatisticalProperty
{
private:
    shared_ptr<System> m_system;
    StatisticalValue<double> m_value;

public:
    PotentialEnergySampler(shared_ptr<System> system);
    virtual void action();
    StatisticalValue<double> value();
};
