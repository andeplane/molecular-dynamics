#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
using std::vector;

#include <potential.h>
#include <topology.h>
#include <systemcell.h>
#include <atom.h>

class Topology;

class System
{
private:
    vector<Potential*> m_potentials;
    Topology m_topology;
    vector<Atom> m_atoms;
    vector<SystemCell> m_cells;
    int m_firstGhostAtomIndex;
    double m_cutoffDistance;
public:
    void initialize(int nodeIndex, int numNodesVector[3], double systemLength[3], double cutoffDistance);
    void resetForces();

    vector<Potential *> &potentials();
    vector<Atom> &atoms();
    void setAtoms(const vector<Atom> &atoms);
    Topology topology() const;
    void setTopology(const Topology &topology);
    vector<SystemCell> cells() const;
    void setCells(const vector<SystemCell> &cells);
    int firstGhostAtomIndex() const;
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    double cutoffDistance() const;
    void setCutoffDistance(double cutoffDistance);
    void updateCells();
    void addPotential(PotentialType type);
};

#endif // SYSTEM_H
