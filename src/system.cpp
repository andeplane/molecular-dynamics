#include "system.h"
#include "atom.h"
#include "systemcell.h"
#include "topology.h"
#include "cmath"
#include "systemcell.h"
#include "potentials/lennardjonespotential.h"
#include <iostream>
using namespace std;


AtomList System::atomList() const
{
    return m_atomList;
}

System::System() :
    m_indexOfNextFreeAtom(0),
    m_firstGhostAtomIndex(-1),
    m_cutoffDistance(1),
    m_isInitialized(false)
{

}

void System::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength, double cutoffDistance)
{
    m_isInitialized = true;
    m_cutoffDistance = cutoffDistance;
    m_firstGhostAtomIndex = -1;
    m_topology.initialize(nodeIndex, numNodesVector, systemLength);

    for(int a=0; a<3; a++) {
        int numCells = ceil(systemLength[a] / cutoffDistance);
        SystemCell::numberOfCellsWithoutGhostCells[a] = numCells;
        SystemCell::numberOfCellsWithGhostCells[a] = numCells+2; // Add two extra ghost cells
        SystemCell::cellLength[a] = systemLength[a] / SystemCell::numberOfCellsWithoutGhostCells[a];
    }

    int numberOfCells = SystemCell::numberOfCellsWithGhostCells[0]*SystemCell::numberOfCellsWithGhostCells[1]*SystemCell::numberOfCellsWithGhostCells[2];
    m_cells.resize(numberOfCells);
}

Atom *System::addAtom()
{
    return m_atomList.addAtom();
}

void System::removeAtom(Atom *atom)
{
    m_atomList.removeAtom(atom);
}

void System::removeAllAtoms()
{
    m_atomList.removeAllAtoms();
}

int System::numberOfAtoms()
{
    return m_atomList.atoms().size();
}

void System::addPotential(PotentialType type) {
    checkIfInitialized();
    if(type == PotentialType::LennardJones) {
        LennardJonesPotential *lennardJones = new LennardJonesPotential();
        potentials().push_back(lennardJones);
    }
}

void System::checkIfInitialized() {
    if(!m_isInitialized) std::cerr << "Error, System object not initialized." << std::endl;
}

void System::updateCells() {
    checkIfInitialized();

    for(SystemCell &cell : m_cells) {
        cell.reset();
    }

    for(Atom *atom : m_atomList.atoms()) {
        int cellIndex = SystemCell::cellIndexForAtom(atom);
        int count = m_cells.size();
        SystemCell &cell = m_cells.at(cellIndex);
        cell.addAtom(atom);
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

vector<Potential*> &System::potentials()
{
    return m_potentials;
}

void System::resetForces() {
    checkIfInitialized();
    for(Atom *atom : m_atomList.atoms()) {
        atom->resetForce();
    }
}

vector<Atom*> &System::atoms()
{
    return m_atomList.atoms();
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
