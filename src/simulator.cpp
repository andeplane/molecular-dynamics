#include "simulator.h"
#include "integrators/integrator.h"
#include "integrators/velocityverlet.h"
#include <iostream>

void Simulator::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength, double timestep, Integrators integrator) {
    m_system.initialize(nodeIndex, numNodesVector, systemLength);
    setIntegrator(integrator);
    m_initialized = true;
    m_timestep = timestep;
}

Simulator::Simulator() :
    m_integrator(NULL),
    m_timestep(0.01),
    m_initialized(false)
{

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

System &Simulator::system()
{
    return m_system;
}

Integrator *Simulator::integrator()
{
    return m_integrator;
}

void Simulator::setIntegrator(Integrators integrator)
{
    if(m_integrator) delete m_integrator;

    if(integrator == Integrators::VelocityVerlet) {
        m_integrator = new VelocityVerlet();
    }
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
    if(!m_initialized) {
        std::cerr << "Simulator not initialized. Remember to call simulator.initialize(...)." << std::endl;
        return;
    }
    m_integrator->integrate(m_system, m_timestep);
    m_timesteps++;          // Increase timestep counter by one
    m_time+= m_timestep;    // Increase time by dt
}
