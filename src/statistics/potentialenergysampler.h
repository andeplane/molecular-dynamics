#pragma once
#include <statistics/statisticalproperty.h>
class System;
namespace Measurements {
    class PotentialEnergySampler : public StatisticalProperty
    {
    private:
        weak_ptr<System> m_system;
        StatisticalValue<double> m_value;

    public:
        static shared_ptr<PotentialEnergySampler> create(shared_ptr<System> system);
        PotentialEnergySampler(shared_ptr<System> system);
        virtual void action();
        StatisticalValue<double> value();
    };
}
