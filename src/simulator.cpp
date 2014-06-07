#include <simulator.h>
#include <integrators/integrator.h>
#include <integrators/velocityverlet.h>
#include <unitconverter.h>

#include <iostream>

void Simulator::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength, double timestep, Integrators integrator) {
    m_system->initialize(nodeIndex, numNodesVector, systemLength);
    setIntegrator(integrator);
    m_initialized = true;
    m_timestep = timestep;
}


vector<shared_ptr<Node> > Simulator::outputs() const
{
    return m_outputs;
}

void Simulator::addOutput(shared_ptr<Node> output)
{
    m_outputs.push_back(output);
}

Simulator::Simulator() :
    m_integrator(NULL),
    m_timestep(UnitConverter::timeFromSI(2e-15)),
    m_initialized(false)
{
    m_system = shared_ptr<System>(new System());
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

shared_ptr<System> Simulator::system()
{
    return m_system;
}

Integrator *Simulator::integrator()
{
    return m_integrator;
}

void Simulator::setIntegrator(Integrators integrator)
{
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

void Simulator::action()
{
    if(!m_initialized) {
        std::cerr << "Simulator not initialized. Remember to call simulator.initialize(...)." << std::endl;
        return;
    }
    m_integrator->integrate(m_system, m_timestep);
    m_timesteps++;          // Increase timestep counter by one
    m_time+= m_timestep;    // Increase time by dt
    for(auto output : m_outputs) { output->step(m_timesteps); }
}
