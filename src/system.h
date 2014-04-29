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

    int m_firstGhostAtomIndex;
    int m_indexOfNextFreeAtom;
    double m_cutoffDistance;
    bool m_isInitialized;
public:
    System();
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength, double cutoffDistance);
    void resetForces();

    vector<Potential *> &potentials();
    vector<Atom*> &atoms();
    Atom *addAtom();
    void removeAtom(Atom *atom);

    Topology topology() const;
    void setTopology(const Topology &topology);
    int firstGhostAtomIndex() const;
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    double cutoffDistance() const;
    void setCutoffDistance(double cutoffDistance);
    void addPotential(PotentialType type);
    AtomManager &atomManager();
    void removeAllAtoms();
    int numberOfAtoms();
protected:
    void checkIfInitialized();
};

#endif // SYSTEM_H
