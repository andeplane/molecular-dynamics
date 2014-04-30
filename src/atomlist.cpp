#include <atomlist.h>
#include <utility>

AtomList::AtomList(int initialAtomCount) :
    m_atomsDirty(false)
{
    m_atoms.reserve(initialAtomCount);
}

AtomList::~AtomList()
{
    m_atoms.clear();
}

int AtomList::numberOfAtoms() const
{
    return m_atoms.size();
}

Atom &AtomList::addAtom(AtomType *atomType)
{
    m_atoms.resize(m_atoms.size()+1); // Default constructor will be called
    m_atoms.back().setType(atomType);
    return m_atoms.back();
}

void AtomList::removeAllAtoms()
{
    m_atoms.clear();
}

void AtomList::iterate(function<void (Atom &atom, const int &atomIndex)> action)
{
    for(int atomIndex=0; atomIndex<m_atoms.size(); atomIndex++) {
        Atom &atom = m_atoms.at(atomIndex);
        action(atom, atomIndex);
    }
}


const vector<Atom> &AtomList::atoms()
{
    return m_atoms;
}

void AtomList::cleanupList() {
    m_atomsDirty = false;
    if(m_atoms.size() == 0) return;
    int numberOfAtoms = m_atoms.size();
    for(int i=0; i<numberOfAtoms; i++) {
        Atom &atom = m_atoms.at(i);
        if(atom.moved()) {
            std::swap(atom,m_atoms.back());
            numberOfAtoms--;
            i--;
        }
    }

    m_atoms.erase(m_atoms.begin()+numberOfAtoms,m_atoms.end());
}
