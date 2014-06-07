#pragma once
#include <memory>
#include <statistics/statisticalproperty.h>
#include <statistics/kineticenergysampler.h>
using std::shared_ptr;

class System;
namespace Measurements {
    class TemperatureSampler : public StatisticalProperty
    {
    private:
        StatisticalValue<double> m_value;
        weak_ptr<Measurements::KineticEnergySampler> m_kineticEnergySampler;
        weak_ptr<System> m_system;
    public:
        static shared_ptr<TemperatureSampler> create(shared_ptr<Measurements::KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system);
        TemperatureSampler(shared_ptr<Measurements::KineticEnergySampler> kineticEnergySampler, shared_ptr<System> system);
        virtual void action();
        StatisticalValue<double> value();
    };
}
