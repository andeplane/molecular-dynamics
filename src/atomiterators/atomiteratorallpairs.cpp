#include <atomiterators/atomiteratorallpairs.h>
#include <cell.h>
#include <atommanager.h>


bool AtomIteratorAllPairs::loopThroughGhosts() const
{
    return m_loopThroughGhosts;
}

void AtomIteratorAllPairs::setLoopThroughGhosts(bool loopThroughGhosts)
{
    m_loopThroughGhosts = loopThroughGhosts;
}

AtomIteratorAllPairs::AtomIteratorAllPairs() :
    m_loopThroughGhosts(false)
{
}

void AtomIteratorAllPairs::iterate(AtomManager &atomManager) {
    atomManager.atoms().iterate([&](Atom &atom1) {
        atomManager.atoms().iterate([&](Atom &atom2) {
            if(atom1.originalUniqueId() <= atom2.originalUniqueId()) return; // Newton's 3rd law
            this->twoParticleAction()(&atom1,&atom2);
        });

        if(m_loopThroughGhosts) {
            atomManager.ghostAtoms().iterate([&](Atom &atom2) {
                if(atom1.originalUniqueId() < atom2.originalUniqueId()) return; // Newton's 3rd law
                this->twoParticleAction()(&atom1,&atom2);
            });
        }
    });
}
