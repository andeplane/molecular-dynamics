#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <system.h>
#include <integrators/integrator.h>
#include <statisticssampler.h>

class Simulator
{
private:
    System m_system;
    Integrator *m_integrator;
    StatisticsSampler m_sampler;
    double m_timestep;
    double m_time;
    int m_timesteps;

public:
    Simulator(int nodeIndex, int numNodesVector[], double systemLength[], double cutoffDistance, double timestep = 0.02);
    System &system();
    Integrator *integrator() const;
    void setIntegrator(Integrator *integrator);
    StatisticsSampler sampler() const;
    void setSampler(const StatisticsSampler &sampler);
    double timestep() const;
    void setTimestep(double timestep);
    double time() const;
    void setTime(double time);
    int timesteps() const;
    void setTimesteps(int timesteps);

    void step();
};

#endif // SIMULATOR_H
