#include <atomlist.h>
#include <utility>
#include <iostream>
using namespace std;

AtomList::AtomList(int initialAtomCount) :
    m_atomsDirty(false),
    m_onAtomMoved(0)
{
    m_atoms.reserve(initialAtomCount);
}

AtomList::~AtomList()
{
    m_atoms.clear();
}

int AtomList::numberOfAtoms()
{
    return atoms().size();
}

Atom &AtomList::addAtom(AtomType *atomType)
{
    m_atoms.push_back(Atom(atomType));
    Atom &atom = m_atoms.back();

    atom.setType(atomType);
    atom.addOnRemoved([&]() {
        // This list should know if an atom has been moved
        m_atomsDirty = true;
    });
    atom.addOnMoved(m_onAtomMoved); // Set default function that will be called on move function on atoms
    return atom;
}

void AtomList::removeAllAtoms()
{
    m_atoms.clear();
}

void AtomList::iterate(function<void (Atom &atom)> action)
{
    for(unsigned long atomIndex=0; atomIndex<m_atoms.size(); atomIndex++) {
        Atom &atom = m_atoms.at(atomIndex);
        action(atom);
    }
}

void AtomList::iterate(function<void (Atom &atom, const int &atomIndex)> action)
{
    for(unsigned long atomIndex=0; atomIndex<m_atoms.size(); atomIndex++) {
        Atom &atom = m_atoms.at(atomIndex);
        action(atom, atomIndex);
    }
}

const vector<Atom> &AtomList::atoms()
{
    if(m_atomsDirty) cleanupList();
    return m_atoms;
}

void AtomList::cleanupList() {
    m_atomsDirty = false;
    if(m_atoms.size() == 0) return;
    int numberOfAtoms = m_atoms.size();

    // Point on the back element so we can switch a moved atom with this one
    vector<Atom>::iterator lastElementIterator = --m_atoms.end();

    for(int atomIndex=0; atomIndex<numberOfAtoms; atomIndex++) {
        Atom &atom = m_atoms.at(atomIndex);
        if(atom.removed()) {
            std::swap(atom,*lastElementIterator);   // Swap this moved element and the back element
            lastElementIterator--;                  // Now, choose the previous element as back element
            numberOfAtoms--;                        // For the for loop range check
            atomIndex--;                            // Re-check this index, a new atom has taken its place
        }
    }

    m_atoms.resize(numberOfAtoms);
}

void AtomList::resetVelocityZero()
{
    iterate([](Atom &atom) {
       atom.setVelocity(0,0,0);
    });
}

void AtomList::resetVelocityMaxwellian(double temperature)
{
    iterate([&](Atom &atom) {
       atom.resetVelocityMaxwellian(temperature);
    });
}

void AtomList::setOnAtomMoved(const function<void ()> &onAtomMoved)
{
    m_onAtomMoved = onAtomMoved;
}

std::ostream& operator<<(std::ostream &stream, AtomList &atomList) {
    stream << "Atom list has " << atomList.numberOfAtoms() << " atoms." << endl;
    atomList.iterate([&](Atom &atom, const int &atomIndex) {
        stream << atomIndex << ": " << atom << endl;
    });

    return stream;
}
