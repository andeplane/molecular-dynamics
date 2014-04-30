#ifndef ATOMLIST_H
#define ATOMLIST_H
#include <vector>
#include <functional>
using std::vector;
using std::function;
#include <atom.h>

class AtomList
{
private:
    vector<Atom> m_atoms;
    bool m_atomsDirty;
    void cleanupList();
public:
    AtomList(int initialAtomCount = 1000);
    ~AtomList();
    int numberOfAtoms() const;
    Atom &addAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    void removeAllAtoms();
    void iterate(function<void(Atom &, const int &)> action);
    vector<Atom> &atoms();
};

#endif // ATOMLIST_H
