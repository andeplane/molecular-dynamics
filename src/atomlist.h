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
    friend std::ostream& operator<<(std::ostream&stream, AtomList&atomList);
    vector<Atom> m_atoms;
    bool m_atomsDirty;
    function<void()> m_onAtomMoved;
    function<void()> m_onAtomRemoved;
public:
    AtomList(int initialAtomCount = 1000);
    ~AtomList();
    int numberOfAtoms();
    Atom &addAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    void removeAllAtoms();
    void iterate(function<void(Atom &, const int &)> action);
    void iterate(function<void(Atom &)> action);
    const vector<Atom> &atoms();
    void setOnAtomMoved(const function<void ()> &onAtomMoved);
    void cleanupList();
    void resetVelocityZero();
    void resetVelocityMaxwellian(double temperature);
};

#endif // ATOMLIST_H
