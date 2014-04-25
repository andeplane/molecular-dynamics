#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
using std::vector;

class SystemCell;
class Atom;
class Potential;
class Topology;

class System
{
private:
    vector<Potential> m_potentials;
    Topology m_topology;
    vector<Atom> m_atoms;
    vector<SystemCell> m_cells;
    int m_firstGhostAtomIndex;
    double m_cutoffDistance;
public:
    System();
    void resetForces();

    Potential potentials() const;
    void setPotential(const Potential &potential);
    vector<Atom> &atoms() const;
    void setAtoms(const vector<Atom> &atoms);
    Topology topology() const;
    void setTopology(const Topology &topology);
    vector<SystemCell> cells() const;
    void setCells(const vector<SystemCell> &cells);
    int firstGhostAtomIndex() const;
    void setFirstGhostAtomIndex(int firstGhostAtomIndex);
    double cutoffDistance() const;
    void setCutoffDistance(double cutoffDistance);
};

#endif // SYSTEM_H
