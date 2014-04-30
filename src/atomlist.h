#ifndef ATOMLIST_H
#define ATOMLIST_H
#include <vector>
using std::vector;

#include <atom.h>

class AtomList
{
private:
    int m_numberOfAtoms;
    vector<Atom> m_atoms;
    bool m_atomsDirty;
    void cleanupList();
public:
    AtomList(int initialAtomCount = 1000);
    ~AtomList();
    int numberOfAtoms() const;
    Atom &addAtom();
    void removeAllAtoms();
};

#endif // ATOMLIST_H
