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
            if(atom1.uniqueId() >= atom2.uniqueId()) return;

            m_twoParticleAction(&atom1,&atom2); // AA
            if(!m_threeParticleAction) return;

            atomManager.atoms().iterate([&](Atom &atom3) {
                if(atom2.uniqueId() >= atom3.uniqueId()) return;
                m_threeParticleAction(&atom1, &atom2, &atom3); // AAA
            });

            if(m_loopThroughGhosts) {
                atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    m_threeParticleAction(&atom1,&atom2, &atom3); // AAG
                });
            }
        });

        if(m_loopThroughGhosts) {
            atomManager.ghostAtoms().iterate([&](Atom &atom2) {
                if(atom1.uniqueId() >= atom2.uniqueId()) return; // Newton's 3rd law
                m_twoParticleAction(&atom1, &atom2); // AG
                if(!m_threeParticleAction) return;

                atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    m_threeParticleAction(&atom1,&atom2, &atom3); // AGG
                });

                atomManager.atoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    m_threeParticleAction(&atom1,&atom2, &atom3); // AGA
                });
            });
        }
    });


}
