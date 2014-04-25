#include "system.h"
#include <atom.h>
#include <potential.h>
#include <systemcell.h>
#include <topology.h>
#include <cmath>
#include <systemcell.h>

void System::initialize(int nodeIndex, int numNodesVector[3], double systemLength[3], double cutoffDistance)
{
    m_cutoffDistance = cutoffDistance;
    m_firstGhostAtomIndex = 0;
    m_topology.initialize(nodeIndex, numNodesVector, systemLength);

    for(int a=0; a<3; a++) {
        int numCells = ceil(systemLength[a] / cutoffDistance);
        SystemCell::numberOfCells[a] = numCells;
        SystemCell::cellLength[a] = systemLength[a] / SystemCell::numberOfCells[a];
    }

    int numberOfCells = SystemCell::numberOfCells[0]*SystemCell::numberOfCells[1]*SystemCell::numberOfCells[2];
    m_cells.resize(numberOfCells);
}

void System::updateCells() {
    for(SystemCell &cell : m_cells) {
        cell.reset();
    }

    for(Atom &atom : m_atoms) {
        SystemCell &cell = m_cells.at(SystemCell::cellIndexForAtom(atom));
        cell.addAtom(&atom);
    }
}

double System::cutoffDistance() const
{
    return m_cutoffDistance;
}

void System::setCutoffDistance(double cutoffDistance)
{
    m_cutoffDistance = cutoffDistance;
}

int System::firstGhostAtomIndex() const
{
    return m_firstGhostAtomIndex;
}

void System::setFirstGhostAtomIndex(int firstGhostAtomIndex)
{
    m_firstGhostAtomIndex = firstGhostAtomIndex;
}

vector<Potential*> System::potentials() const
{
    return m_potentials;
}

void System::resetForces() {
    for(Atom atom : m_atoms) {
        atom.resetForce();
    }
}

vector<Atom> &System::atoms()
{
    return m_atoms;
}

void System::setAtoms(const vector<Atom> &atoms)
{
    m_atoms = atoms;
}

Topology System::topology() const
{
    return m_topology;
}

void System::setTopology(const Topology &topology)
{
    m_topology = topology;
}

vector<SystemCell> System::cells() const
{
    return m_cells;
}

void System::setCells(const vector<SystemCell> &cells)
{
    m_cells = cells;
}
