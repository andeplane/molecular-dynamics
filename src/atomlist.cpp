#include <atomlist.h>
#include <atom.h>
#include <iostream>
using std::cerr; using std::endl;

AtomList::AtomList(int initialAtomCount) :
    m_nextFreeIndex(0)
{
    m_allAtoms.resize(initialAtomCount);
}

vector<Atom *> &AtomList::atoms()
{
    return m_atoms;
}

int AtomList::indexInAllAtoms(Atom *atom)
{
    pair<int,int> &indexPair = m_indexMap.find(atom)->second;
    return indexPair.first;
}

void AtomList::increaseNumberOfAtoms()
{
    int newSize = 2*(m_allAtoms.size()+1);

    vector<int> indicesOfAtomsInUse;
    indicesOfAtomsInUse.reserve(m_atoms.size());
    for(Atom *atom : m_atoms) {
        indicesOfAtomsInUse.push_back(indexInAllAtoms(atom));
    }

    m_allAtoms.resize(newSize);
    m_atoms.clear();
    m_indexMap.clear();

    for(int indexInAllAtoms : indicesOfAtomsInUse) {
        m_indexMap.insert(pair<Atom*,pair<int,int>>(&m_allAtoms.at(indexInAllAtoms), pair<int,int>(indexInAllAtoms, m_atoms.size())));
        m_atoms.push_back(&m_allAtoms.at(indexInAllAtoms));
    }

    indicesOfAtomsInUse.clear();
}

int AtomList::indexInAtoms(Atom *atom)
{
    pair<int,int> &indexPair = m_indexMap.find(atom)->second;
    return indexPair.second;
}

Atom *AtomList::addAtom()
{
    if(m_nextFreeIndex>=m_allAtoms.size()) {
        increaseNumberOfAtoms();
    }

    int indexInAtoms = m_atoms.size();
    int indexInAllAtoms = m_nextFreeIndex;

    Atom *nextAtom = &m_allAtoms.at(indexInAllAtoms);

    m_indexMap.insert(pair<Atom*,pair<int,int>>(nextAtom, pair<int,int>(indexInAllAtoms, indexInAtoms)));
    m_atoms.push_back(nextAtom);
    m_nextFreeIndex++;
    return nextAtom;
}

void AtomList::removeAtom(Atom *atom)
{
    if(m_indexMap.count(atom) == 0) {
        cerr << "Tried to remove atom not in used. Check your code." << endl;
        return;
    }

    int removedAtomIndexAtoms = indexInAtoms(atom);
    int removedAtomIndexAllAtoms = indexInAllAtoms(atom);
    int lastUsedAtomIndexAllAtoms = --m_nextFreeIndex;

    std::swap(m_allAtoms.at(removedAtomIndexAllAtoms), m_allAtoms.at(lastUsedAtomIndexAllAtoms));
    m_atoms.pop_back();

    // allAtoms.at(removedAtomIndexAllAtoms) = allAtoms.at(lastUsedAtomIndexAllAtoms);
    Atom *removedAtomAllAtoms = &m_allAtoms.at(removedAtomIndexAllAtoms);
    m_indexMap.insert(pair<Atom*,pair<int,int>>(removedAtomAllAtoms, pair<int,int>(removedAtomIndexAllAtoms, removedAtomIndexAtoms)));

    Atom *lastUsedAtomAllAtoms = &m_allAtoms.at(lastUsedAtomIndexAllAtoms);
    m_indexMap.erase(m_indexMap.find(lastUsedAtomAllAtoms)); // Remove the pointer of the atom that is not used anymore
}

int AtomList::atomCapacity()
{
    return m_allAtoms.size();
}

void AtomList::removeAllAtoms()
{
    m_indexMap.clear();
    m_atoms.clear();
    m_nextFreeIndex = 0;
}
