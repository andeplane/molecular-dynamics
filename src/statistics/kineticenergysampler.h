#pragma once
#include <statistics/statisticalproperty.h>

class KineticEnergySampler : public StatisticalProperty
{
private:
    StatisticalValue<double> m_value;
    System *m_system;
public:
    KineticEnergySampler(System *system);
    StatisticalValue<double> value();
    virtual void action();
};
