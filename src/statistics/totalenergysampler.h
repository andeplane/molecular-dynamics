#pragma once
#include <statistics/statisticalproperty.h>
class KineticEnergySampler;
class PotentialEnergySampler;

class TotalEnergySampler : public StatisticalProperty
{
private:
    StatisticalValue<double> m_value;
    shared_ptr<KineticEnergySampler> m_kineticEnergySampler;
    shared_ptr<PotentialEnergySampler> m_potentialEnergySampler;

public:
    TotalEnergySampler(shared_ptr<KineticEnergySampler> kineticEnergySampler, shared_ptr<PotentialEnergySampler> potentialEnergySampler);
    virtual void action();
    StatisticalValue<double> value();
};
