#ifndef SIMULATOR_H
#define SIMULATOR_H
class System;
class Integrator;
class StatisticsSampler;

class Simulator
{
private:
    System m_system;
    Integrator m_integrator;
    StatisticsSampler m_sampler;
    double m_timestep;
    double m_time;
    int m_timesteps;

public:
    Simulator(double timestep = 0.02);
    System system() const;
    void setSystem(const System &system);
    Integrator integrator() const;
    void setIntegrator(const Integrator &integrator);
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
