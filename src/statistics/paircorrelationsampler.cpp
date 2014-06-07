#include <statistics/paircorrelationsampler.h>
#include <system.h>
#include <statistics/neighborlist.h>
PairCorrelationSampler::PairCorrelationSampler(shared_ptr<AtomType> atomType1, shared_ptr<AtomType> atomType2, shared_ptr<System> system, double maxDistance) :
    m_atomType1(atomType1),
    m_atomType2(atomType2),
    m_system(system)
{
    m_neighborList = shared_ptr<NeighborList>(new NeighborList(system, maxDistance));
}

void PairCorrelationSampler::action()
{

}
