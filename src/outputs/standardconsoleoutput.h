#pragma once
#include <node.h>
#include <statistics/statistics.h>
#include <system.h>

class StandardConsoleOutput : public Node
{
private:
    shared_ptr<Measurements::TotalEnergySampler> m_totalEnergySampler;
    shared_ptr<Measurements::TemperatureSampler> m_temperatureSampler;
    shared_ptr<System> m_system;
    int m_frequency;

public:
    static shared_ptr<StandardConsoleOutput> create(shared_ptr<Measurements::TotalEnergySampler> totalEnergySampler, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system, int frequency = 100);
    StandardConsoleOutput(shared_ptr<Measurements::TotalEnergySampler> totalEnergySampler, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system, int frequency = 100);
    virtual void action();
    int frequency() const;
    void setFrequency(int frequency);
};
