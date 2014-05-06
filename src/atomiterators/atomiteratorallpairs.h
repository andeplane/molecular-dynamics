#pragma once
#include <atomiterators/atomiterator.h>

class AtomIteratorAllPairs : public AtomIterator
{
private:
    bool m_loopThroughGhosts;
public:
    AtomIteratorAllPairs();
    virtual void iterate(AtomManager &atomManager);
    bool loopThroughGhosts() const;
    void setLoopThroughGhosts(bool loopThroughGhosts);
};
