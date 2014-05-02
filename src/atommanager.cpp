#include <atommanager.h>
#include <atom.h>
#include <iostream>
#include <cmath>
#include <topology.h>

using std::cerr; using std::endl; using std::cout;

AtomManager::~AtomManager()
{
    removeAllAtoms();
    m_cellData.cells.clear();
    m_mpiReceiveBuffer.clear();
    m_mpiSendBuffer.clear();
}

AtomManager::AtomManager() :
    m_cellStructureDirty(true),
    m_cellDataDirty(true),
    m_ghostAtomsDirty(true),
    m_topology(0),
    m_updatingGhostAtoms(false),
    m_ghostAtomsEnabled(true)
{
    m_cellData.cutoffDistance = INFINITY;
    m_cellData.initialized = false;

    // Default behaviour when atoms move
    m_atoms.setOnAtomMoved([&]() {
        m_cellDataDirty = true;
        m_ghostAtomsDirty = true;
    });
}

CellData &AtomManager::cellData()
{
    if(m_cellDataDirty) updateCellList();
    return m_cellData;
}

int AtomManager::numberOfAtoms()
{
    return atoms().numberOfAtoms();
}

int AtomManager::numberOfGhostAtoms()
{
    return ghostAtoms().numberOfAtoms();
}

Atom &AtomManager::addAtom(AtomType *atomType)
{
    return m_atoms.addAtom(atomType);
}

Atom &AtomManager::addGhostAtom(AtomType *atomType)
{
    Atom &atom = m_ghostAtoms.addAtom(atomType);
    atom.setGhost(true);
    return atom;
}

void AtomManager::removeGhostAtoms()
{
    m_ghostAtoms.removeAllAtoms();
    m_ghostAtomsDirty = true;
}

void AtomManager::removeAllAtoms()
{
    m_atoms.removeAllAtoms();
    m_ghostAtoms.removeAllAtoms();
    m_ghostAtomsDirty = true;
}

void AtomManager::setCutoffDistance(double cutoffDistance) {
    if(m_cellData.cutoffDistance != cutoffDistance) {
        m_ghostAtomsDirty = true;
        m_cellStructureDirty = true;
        m_cellData.cutoffDistance = cutoffDistance;
    }
}

double AtomManager::cutoffDistance()
{
    return m_cellData.cutoffDistance;
}

void AtomManager::setSystemLength(vector<double> &systemLength) {
    if(systemLength.at(0) != m_cellData.systemLength[0] || systemLength.at(1) != m_cellData.systemLength[1] || systemLength.at(2) != m_cellData.systemLength[2]) {
        m_ghostAtomsDirty = true;
        m_cellStructureDirty = true;
        m_cellData.systemLength[0] = systemLength.at(0);
        m_cellData.systemLength[1] = systemLength.at(1);
        m_cellData.systemLength[2] = systemLength.at(2);
    }

    m_cellData.initialized = true;
}


AtomList &AtomManager::atoms()
{
    return m_atoms;
}

void AtomManager::updateGhostAtoms() {
    m_updatingGhostAtoms = true;
    m_topology->copyGhostAtomsWithMPI(*this);
    m_updatingGhostAtoms = false;
    m_ghostAtomsDirty = false;
}

AtomList &AtomManager::ghostAtoms()
{
    if(m_ghostAtomsEnabled) {
        if(m_ghostAtomsDirty && !m_updatingGhostAtoms) updateGhostAtoms();
    } else removeGhostAtoms();

    return m_ghostAtoms;
}

void AtomManager::setTopology(Topology *topology)
{
    m_topology = topology;
}

bool AtomManager::ghostAtomsEnabled() const
{
    return m_ghostAtomsEnabled;
}

void AtomManager::setGhostAtomsEnabled(bool ghostAtomsEnabled)
{
    m_ghostAtomsEnabled = ghostAtomsEnabled;
}
void AtomManager::updateCellStructure() {
    if(!m_cellData.initialized) {
        std::cerr << "atomManager.cellData not initialized. We don't know the systemLength." << std::endl;
        return;
    }

    if(isinf(m_cellData.cutoffDistance)) {
        m_cellData.numberOfCellsWithoutGhostCells[0] = 1;
        m_cellData.numberOfCellsWithoutGhostCells[1] = 1;
        m_cellData.numberOfCellsWithoutGhostCells[2] = 1;
    } else {
        m_cellData.numberOfCellsWithoutGhostCells[0] = ceil( m_cellData.systemLength[0] / m_cellData.cutoffDistance);
        m_cellData.numberOfCellsWithoutGhostCells[1] = ceil( m_cellData.systemLength[1] / m_cellData.cutoffDistance);
        m_cellData.numberOfCellsWithoutGhostCells[2] = ceil( m_cellData.systemLength[2] / m_cellData.cutoffDistance);
    }

    m_cellData.numberOfCellsWithGhostCells[0] = m_cellData.numberOfCellsWithoutGhostCells[0]+2;
    m_cellData.numberOfCellsWithGhostCells[1] = m_cellData.numberOfCellsWithoutGhostCells[1]+2;
    m_cellData.numberOfCellsWithGhostCells[2] = m_cellData.numberOfCellsWithoutGhostCells[2]+2;

    unsigned long numCellsTotal = m_cellData.numberOfCellsWithGhostCells[0]*m_cellData.numberOfCellsWithGhostCells[1]*m_cellData.numberOfCellsWithGhostCells[2];

    if(m_cellData.cells.size() < numCellsTotal) {
        m_cellData.cells.resize(numCellsTotal);
    }

    for(unsigned long cellIndex=0; cellIndex<m_cellData.cells.size(); cellIndex++) {
        Cell &cell = m_cellData.cells.at(cellIndex);
        cell.setCellIndex(cellIndex);
    }

    m_cellData.cellLength[0] = m_cellData.systemLength[0] / m_cellData.numberOfCellsWithoutGhostCells[0];
    m_cellData.cellLength[1] = m_cellData.systemLength[1] / m_cellData.numberOfCellsWithoutGhostCells[1];
    m_cellData.cellLength[2] = m_cellData.systemLength[2] / m_cellData.numberOfCellsWithoutGhostCells[2];

    m_cellStructureDirty = false;
    m_cellDataDirty = true;
}

void AtomManager::updateCellList() {
    if(m_cellStructureDirty) updateCellStructure();
    m_cellDataDirty = false;

    for(Cell &cell : m_cellData.cells) {
        cell.reset();
    }

    CellData &cellData = m_cellData;
    atoms().iterate([&](Atom &atom) {
        int cellIndex = Cell::cellIndexForAtom(atom, cellData);
        Cell &cell = cellData.cells.at(cellIndex);
        cell.addAtom(&atom);
    });

    ghostAtoms().iterate([&](Atom &atom) {
        int cellIndex = Cell::cellIndexForAtom(atom, cellData);
        Cell &cell = cellData.cells.at(cellIndex);
        cell.addAtom(&atom);
    });

}

std::ostream& operator<<(std::ostream &stream, AtomManager &atomManager) {
    return stream << "Atom manager: " << std::endl << "Ghost atoms: " << atomManager.ghostAtoms() << endl << "Atoms: " << atomManager.atoms();
}
