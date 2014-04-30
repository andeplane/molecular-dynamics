#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <topology.h>
#include <atom.h>
#include <atommanager.h>

class Topology;

class System
{
private:
    vector<Potential*> m_potentials;
    Topology m_topology;
    AtomManager m_atomManager;

    int m_numberOfGhostAtoms;
    bool m_isInitialized;
public:
    System();
    ~System();
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength);
    void resetForces();

    vector<Potential *> &potentials();
    vector<Atom*> &atoms();
    Atom *addAtom();
    Atom *addGhostAtom();
    void removeAtom(Atom *atom);

    Topology &topology();
    AtomManager &atomManager();
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    void addPotential(PotentialType type);
    void removeAllAtoms();
    int numberOfAtoms();
    int numberOfGhostAtoms();
    void setSystemLength(vector<double> &systemLength);
protected:
    void checkIfInitialized();
};

#endif // SYSTEM_H
