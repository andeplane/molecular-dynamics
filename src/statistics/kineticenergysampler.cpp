#include <statistics/kineticenergysampler.h>
#include <atom.h>
#include <atommanager.h>

StatisticalValue<double> KineticEnergySampler::value()
{
    return m_value;
}

void KineticEnergySampler::action()
{
    double kineticEnergy = 0;
    m_system->atomManager().atoms().iterate([&] (Atom &atom) {
        kineticEnergy += 0.5*atom.type()->mass()*(atom.velocity[0]*atom.velocity[0] + atom.velocity[1]*atom.velocity[1] + atom.velocity[2]*atom.velocity[2]);
    });
    m_value.addValue(kineticEnergy);
}

KineticEnergySampler::KineticEnergySampler(shared_ptr<System> system) :
    m_system(system)
{

}
