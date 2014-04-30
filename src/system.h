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

    int m_indexOfFirstGhostAtom;
    bool m_isInitialized;
public:
    System();
    ~System();
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength);
    void resetForces();

    vector<Potential *> &potentials();
    vector<Atom*> &atoms();
    Atom *addAtom();
    void removeAtom(Atom *atom);

    Topology &topology();
    int firstGhostAtomIndex() const;
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    void addPotential(PotentialType type);
    AtomManager &atomManager();
    void removeAllAtoms();
    int numberOfAtoms();
    int numberOfGhostAtoms();
    void setSystemLength(vector<double> &systemLength);

protected:
    void checkIfInitialized();
};

#endif // SYSTEM_H
