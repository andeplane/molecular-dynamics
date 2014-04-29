#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
using std::vector;

#include <potentials/potential.h>
#include <topology.h>
#include <systemcell.h>
#include <atom.h>
#include <atomlist.h>

class Topology;

class System
{
private:
    vector<Potential*> m_potentials;
    Topology m_topology;
    AtomList m_atomList;
    vector<SystemCell> m_cells;
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
    vector<SystemCell> cells() const;
    int firstGhostAtomIndex() const;
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    double cutoffDistance() const;
    void setCutoffDistance(double cutoffDistance);
    void updateCells();
    void addPotential(PotentialType type);
    AtomList atomList() const;
    void removeAllAtoms();
    int numberOfAtoms();
protected:
    void checkIfInitialized();
};

#endif // SYSTEM_H
