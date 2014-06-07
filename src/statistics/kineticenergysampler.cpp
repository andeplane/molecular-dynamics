#include <statistics/kineticenergysampler.h>
#include <atom.h>
#include <atommanager.h>
#include <system.h>
#include <includes.h>
using namespace Measurements;
StatisticalValue<double> KineticEnergySampler::value()
{
    return m_value;
}

void KineticEnergySampler::action()
{
    double kineticEnergy = 0;
    m_system.lock()->atomManager().atoms().iterate([&] (Atom &atom) {
        kineticEnergy += 0.5*atom.type()->mass()*(atom.velocity[0]*atom.velocity[0] + atom.velocity[1]*atom.velocity[1] + atom.velocity[2]*atom.velocity[2]);
    });
    m_value.addValue(kineticEnergy);
}

shared_ptr<KineticEnergySampler> KineticEnergySampler::create(shared_ptr<System> system)
{
    return make_shared<KineticEnergySampler>(system);
}

KineticEnergySampler::KineticEnergySampler(shared_ptr<System> system) :
    m_system(system)
{

}
