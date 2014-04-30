#ifndef ATOMITERATORDEFAULT_H
#define ATOMITERATORDEFAULT_H
#include <atomiterators/atomiterator.h>

class AtomIteratorDefault : public AtomIterator
{
public:
    AtomIteratorDefault();
    virtual void iterate(AtomManager &atomManager);
};

#endif // ATOMITERATORDEFAULT_H
