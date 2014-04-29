#ifndef ATOMLIST_H
#define ATOMLIST_H

#include <vector>
#include <map>
using std::vector; using std::map; using std::pair;

class Atom;

class AtomList
{
private:
    vector<Atom> m_allAtoms;    // All atoms
    vector<Atom*> m_atoms;      // Atoms in use
    int m_nextFreeIndex;
    map<Atom*,pair<int,int> > m_indexMap;
    int indexInAtoms(Atom *atom);
    int indexInAllAtoms(Atom *atom);
    void increaseNumberOfAtoms();

public:
    AtomList(int initialAtomCount = 1000);
    vector<Atom *> &atoms();
    Atom *addAtom();
    void removeAtom(Atom *atom);
    int atomCapacity();
    void removeAllAtoms();
};

#endif // ATOMLIST_H
