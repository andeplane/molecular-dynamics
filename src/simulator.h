#pragma once
#include <system.h>
#include <integrators/integrator.h>
#include <statisticssampler.h>
#include <node.h>

class Simulator : public Node
{
private:
    shared_ptr<System> m_system;
    Integrator *m_integrator;
    vector<shared_ptr<Node>> m_outputs;
    double m_timestep;
    double m_time;
    int m_timesteps;
    bool m_initialized;

public:
    Simulator();
    shared_ptr<System> system();
    void initialize(int nodeIndex, vector<int> numNodesVector, vec3 systemLength, double timestep = 0.01, Integrators integrator = Integrators::VelocityVerlet);
    Integrator *integrator();
    void setIntegrator(Integrators integrator);
    double timestep() const;
    void setTimestep(double timestep);
    double time() const;
    void setTime(double time);
    int timesteps() const;
    void setTimesteps(int timesteps);
    virtual void action();
    vector<shared_ptr<Node> > outputs() const;
    void addOutput(shared_ptr<Node> output);
};
