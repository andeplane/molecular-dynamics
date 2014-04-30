#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <topology.h>
#include <atom.h>
#include <atommanager.h>

class System
{
private:
    vector<Potential*> m_potentials;
    Topology m_topology;
    AtomManager m_atomManager;
    bool m_isInitialized;
public:
    System();
    ~System();
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength);

    vector<Potential *> &potentials();
    Atom &addAtom();
    Atom &addGhostAtom();
    void removeAllAtoms();
    void removeGhostAtoms();
    void setSystemLength(vector<double> &systemLength);

    Topology &topology();
    AtomManager &atomManager();
    void addPotential(PotentialType type);
    int numberOfAtoms();
    int numberOfGhostAtoms();
protected:
    void checkIfInitialized();
};

#endif // SYSTEM_H
