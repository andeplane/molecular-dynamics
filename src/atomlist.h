#pragma once
#include <vector>
#include <functional>
#include <unordered_map>
using std::vector; using std::function; using std::unordered_map; using std::pair;

#include <particles/atom.h>

class AtomList
{
private:
    vector<Atom> m_atoms;
    bool m_atomsDirty;
    unordered_map<unsigned long, unsigned long> m_indexMap;
    bool m_isIterating;
    friend std::ostream& operator<<(std::ostream&stream, AtomList&atomList);
    function<void()> m_onAtomMoved;
    function<void()> m_onAtomRemoved;
    vector<Atom> m_tempAtoms; // If we add atoms while looping
    void rebuildIndexMap();
    void moveTempAtomsToList();
public:
    AtomList(int initialAtomCount = 1000);
    ~AtomList();
    int numberOfAtoms();
    Atom &addAtom(shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
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
