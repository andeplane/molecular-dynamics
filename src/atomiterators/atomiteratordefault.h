#ifndef ATOMITERATORDEFAULT_H
#define ATOMITERATORDEFAULT_H
#include <atomiterators/atomiterator.h>

class AtomIteratorDefault : public AtomIterator
{
private:
    double m_maximumNeighborDistance;
public:
    AtomIteratorDefault();
    virtual void iterate(AtomManager &atomManager);
    double maximumNeighborDistance() const;
    void setMaximumNeighborDistance(double maximumNeighborDistance);
};

#endif // ATOMITERATORDEFAULT_H
