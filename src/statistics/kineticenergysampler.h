#pragma once
#include <statistics/statisticalproperty.h>
namespace Measurements {
    class KineticEnergySampler : public StatisticalProperty
    {
    private:
        StatisticalValue<double> m_value;
        weak_ptr<System> m_system;
    public:
        static shared_ptr<KineticEnergySampler> create(shared_ptr<System> system);
        KineticEnergySampler(shared_ptr<System> system);
        StatisticalValue<double> value();
        virtual void action();
    };
}
