#pragma once
#include <statistics/statisticalproperty.h>
#include <atomtype.h>

class System; class NeighborList;
class PairCorrelationSampler : public StatisticalProperty
{
private:
    shared_ptr<AtomType> m_atomType1;
    shared_ptr<AtomType> m_atomType2;
    weak_ptr<System> m_system;
    StatisticalValue<int> m_value;
    int m_numberOfBins;
    shared_ptr<NeighborList> m_neighborList;
public:
    PairCorrelationSampler(shared_ptr<AtomType> atomType1, shared_ptr<AtomType> atomType2, shared_ptr<System> system, double maxDistance, int numberOfBins);
    virtual void action();
    double maxDistance();
    int numberOfBins() const;
    void setNumberOfBins(int numberOfBins);
};
