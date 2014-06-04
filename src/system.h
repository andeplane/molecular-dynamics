#pragma once
#include <vector>
#include <node.h>

using std::vector;

#include <potentials/potentials.h>
#include <topology.h>
#include <atom.h>
#include <atommanager.h>

class System : public Node
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
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength);

    vector<Potential *> &potentials();
    Atom &addAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom), vector<double> position = {0,0,0});
    Atom &addGhostAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    void removeAllAtoms();
    void removeGhostAtoms();
    void setSystemLength(vector<double> systemLength);

    Topology &topology();
    AtomManager &atomManager();
    Potential *addPotential(PotentialType type);
    int numberOfAtoms();
    int numberOfGhostAtoms();
    vector<double> systemLength();
    void removeTotalMomentum();
protected:
    void checkIfInitialized();
};
