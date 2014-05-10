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
    // cout << atomManager << endl;
    atomManager.atoms().iterate([&](Atom &atom1) {
        atomManager.atoms().iterate([&](Atom &atom2) {
            if(atom1.uniqueId() >= atom2.uniqueId()) return;

            m_twoParticleAction(&atom1,&atom2);
            if(!m_threeParticleAction) return;

            atomManager.atoms().iterate([&](Atom &atom3) {
                if(atom2.uniqueId() >= atom3.uniqueId()) return;
                // cout << "Triplet: (" << atom1.originalUniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.originalUniqueId() << (atom2.ghost() ? "*" : "") << "," << atom3.originalUniqueId() << (atom3.ghost() ? "*" : "") << ")  --- unique ids (" << atom1.uniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.uniqueId() << (atom2.ghost() ? "*" : "") << ", " << atom3.uniqueId() << (atom3.ghost() ? "*" : "") << ")" << endl;
                m_threeParticleAction(&atom1, &atom2, &atom3);
            });

            if(m_loopThroughGhosts) {
                atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    // cout << "Triplet: (" << atom1.originalUniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.originalUniqueId() << (atom2.ghost() ? "*" : "") << "," << atom3.originalUniqueId() << (atom3.ghost() ? "*" : "") << ")  --- unique ids (" << atom1.uniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.uniqueId() << (atom2.ghost() ? "*" : "") << ", " << atom3.uniqueId() << (atom3.ghost() ? "*" : "") << ")" << endl;
                    m_threeParticleAction(&atom1,&atom2, &atom3);
                });
            }
        });

        if(m_loopThroughGhosts) {
            atomManager.ghostAtoms().iterate([&](Atom &atom2) {
                if(atom1.uniqueId() >= atom2.uniqueId()) return; // Newton's 3rd law
                m_twoParticleAction(&atom1, &atom2);
                if(!m_threeParticleAction) return;

                atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    // cout << "Triplet: (" << atom1.originalUniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.originalUniqueId() << (atom2.ghost() ? "*" : "") << "," << atom3.originalUniqueId() << (atom3.ghost() ? "*" : "") << ")  --- unique ids (" << atom1.uniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.uniqueId() << (atom2.ghost() ? "*" : "") << ", " << atom3.uniqueId() << (atom3.ghost() ? "*" : "") << ")" << endl;
                    m_threeParticleAction(&atom1,&atom2, &atom3);
                });

                atomManager.atoms().iterate([&](Atom &atom3) {
                    if(atom2.uniqueId() >= atom3.uniqueId()) return;
                    // cout << "Triplet: (" << atom1.originalUniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.originalUniqueId() << (atom2.ghost() ? "*" : "") << "," << atom3.originalUniqueId() << (atom3.ghost() ? "*" : "") << ")  --- unique ids (" << atom1.uniqueId() << (atom1.ghost() ? "*" : "") << "," << atom2.uniqueId() << (atom2.ghost() ? "*" : "") << ", " << atom3.uniqueId() << (atom3.ghost() ? "*" : "") << ")" << endl;
                    m_threeParticleAction(&atom1,&atom2, &atom3);
                });
            });
        }
    });


}
