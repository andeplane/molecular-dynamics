#pragma once
#include <vector>
#include <functional>
#include <unordered_map>
using std::vector; using std::function; using std::unordered_map; using std::pair;

#include <atom.h>

class AtomList
{
private:
    friend std::ostream& operator<<(std::ostream&stream, AtomList&atomList);
    vector<Atom> m_atoms;
    bool m_atomsDirty;
    function<void()> m_onAtomMoved;
    function<void()> m_onAtomRemoved;
    unordered_map<unsigned long, unsigned long> m_indexMap;
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
    bool containsAtomWithUniqueId(unsigned long uniqueId);
    Atom &getAtomByUniqueId(unsigned long uniqueId);
};
