#pragma once
#include <atomiterators/atomiterator.h>

class AtomIteratorDefault : public AtomIterator
{
private:
    double m_maximumNeighborDistance;
    bool m_createNeighborList;
public:
    AtomIteratorDefault();
    virtual void iterate(AtomManager &atomManager);
    double maximumNeighborDistance() const;
    void setMaximumNeighborDistance(double maximumNeighborDistance);
    bool createNeighborList() const;
    void setCreateNeighborList(bool createNeighborList);
    virtual void setThreeParticleAction(const function<void (Atom *atom1, Atom *atom2, Atom *atom3)> &threeParticleAction);
};
