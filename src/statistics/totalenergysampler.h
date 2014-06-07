#pragma once
#include <statistics/statisticalproperty.h>
#include <statistics/kineticenergysampler.h>
#include <statistics/potentialenergysampler.h>

namespace Measurements {
    class TotalEnergySampler : public StatisticalProperty
    {
    private:
        StatisticalValue<double> m_value;
        weak_ptr<Measurements::KineticEnergySampler> m_kineticEnergySampler;
        weak_ptr<Measurements::PotentialEnergySampler> m_potentialEnergySampler;

    public:
        static shared_ptr<TotalEnergySampler> create(shared_ptr<Measurements::KineticEnergySampler> kineticEnergySampler, shared_ptr<Measurements::PotentialEnergySampler> potentialEnergySampler);
        TotalEnergySampler(shared_ptr<Measurements::KineticEnergySampler> kineticEnergySampler, shared_ptr<Measurements::PotentialEnergySampler> potentialEnergySampler);
        virtual void action();
        StatisticalValue<double> value();
    };
}
