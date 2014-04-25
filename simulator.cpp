#include "simulator.h"
#include <integrator.h>

Simulator::Simulator(int nodeIndex, int numNodesVector[3], double systemLength[3], double cutoffDistance, double timestep)
{
    m_system.initialize(nodeIndex, numNodesVector, systemLength, cutoffDistance);
}

double Simulator::timestep() const
{
    return m_timestep;
}

void Simulator::setTimestep(double timestep)
{
    m_timestep = timestep;
}

double Simulator::time() const
{
    return m_time;
}

void Simulator::setTime(double time)
{
    m_time = time;
}

int Simulator::timesteps() const
{
    return m_timesteps;
}

void Simulator::setTimesteps(int timesteps)
{
    m_timesteps = timesteps;
}

System Simulator::system() const
{
    return m_system;
}

void Simulator::setSystem(const System &system)
{
    m_system = system;
}

Integrator *Simulator::integrator() const
{
    return m_integrator;
}

void Simulator::setIntegrator(Integrator *integrator)
{
    m_integrator = integrator;
}

StatisticsSampler Simulator::sampler() const
{
    return m_sampler;
}

void Simulator::setSampler(const StatisticsSampler &sampler)
{
    m_sampler = sampler;
}

void Simulator::step() {
    m_integrator->integrate(m_system, m_timestep);
    m_timesteps++;          // Increase timestep counter by one
    m_time+= m_timestep;    // Increase time by dt
}
