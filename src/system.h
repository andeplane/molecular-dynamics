#pragma once
#include <vector>
#include <node.h>

using std::vector;

#include <potentials/potentials.h>
#include <topology.h>
#include <particles/atom.h>
#include <atommanager.h>
#include <utils/vec3.h>

class System : public Node, public std::enable_shared_from_this<System>
{
private:
    friend std::ostream& operator<<(std::ostream&stream, System &system);
    vector<Potential*> m_potentials;
    Topology m_topology;
    AtomManager m_atomManager;
    bool m_isInitialized;
public:
    System();
    ~System();
    void initialize(int nodeIndex, vector<int> numNodesVector, CompPhys::vec3 systemLength);

    vector<Potential *> &potentials();
    Atom &addAtom(shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom), vector<double> position = {0,0,0});
    Atom &addGhostAtom(shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    void removeAllAtoms();
    void removeGhostAtoms();
    void setSystemLength(CompPhys::vec3 systemLength);

    Topology &topology();
    AtomManager &atomManager();
    Potential *addPotential(PotentialType type);
    int numberOfAtoms();
    int numberOfGhostAtoms();
    CompPhys::vec3 systemLength() const;
    void removeTotalMomentum();
protected:
    void checkIfInitialized();
};
