#include <atomlist.h>
#include <utility>

AtomList::AtomList(int initialAtomCount) :
    m_numberOfAtoms(0),
    m_atomsDirty(false)
{
    m_atoms.resize(initialAtomCount);
}

AtomList::~AtomList()
{
    m_atoms.clear();
}

int AtomList::numberOfAtoms() const
{
    return m_numberOfAtoms;
}

Atom &AtomList::addAtom()
{
    if(m_numberOfAtoms>=m_atoms.size()) {
        int newNumberOfAtoms = 2*(m_atoms.size()+1);
        m_atoms.resize(newNumberOfAtoms);
    }

    Atom &atom = m_atoms.at(m_numberOfAtoms++);
    return atom;
}

void AtomList::removeAllAtoms()
{
    m_numberOfAtoms = 0;
}

void AtomList::cleanupList() {
    m_atomsDirty = false;
    if(m_numberOfAtoms == 0) return;
    Atom &lastAtom = m_atoms.at(m_numberOfAtoms-1);

    for(int i=0; i<m_numberOfAtoms; i++) {
        Atom &atom = m_atoms.at(i);
        if(atom.moved()) {
            std::swap(atom,lastAtom);
            m_numberOfAtoms--;
            i--;
        }
    }
}
