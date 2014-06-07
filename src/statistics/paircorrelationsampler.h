#pragma once
#include <statistics/statisticalproperty.h>
#include <atomtype.h>

class System; class NeighborList;
class PairCorrelationSampler : public StatisticalProperty
{
private:
    shared_ptr<AtomType> m_atomType1;
    shared_ptr<AtomType> m_atomType2;
    shared_ptr<System> m_system;
    StatisticalValue<double> m_value;
    shared_ptr<NeighborList> m_neighborList;
public:
    PairCorrelationSampler(shared_ptr<AtomType> atomType1, shared_ptr<AtomType> atomType2, shared_ptr<System> system, double maxDistance);
    virtual void action();
};
