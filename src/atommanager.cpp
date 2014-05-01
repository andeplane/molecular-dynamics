#include <atommanager.h>
#include <atom.h>
#include <iostream>
#include <cmath>

using std::cerr; using std::endl; using std::cout;

AtomManager::~AtomManager()
{

}

AtomManager::AtomManager() :
    m_cellStructureDirty(true),
    m_cellDataDirty(true)
{
    m_cellData.cutoffDistance = INFINITY;
    m_cellData.initialized = false;

    // Default behaviour when atoms move
    m_atoms.setOnAtomMoved([&]() {
        m_cellDataDirty = true;
    });
}

CellData &AtomManager::cellData()
{
    if(m_cellDataDirty) updateCellList();
    return m_cellData;
}

int AtomManager::numberOfAtoms()
{
    return m_atoms.numberOfAtoms();
}

int AtomManager::numberOfGhostAtoms()
{
    return m_ghostAtoms.numberOfAtoms();
}

Atom &AtomManager::addAtom()
{
    return m_atoms.addAtom();
}

Atom &AtomManager::addGhostAtom()
{
    Atom &atom = m_ghostAtoms.addAtom();
    atom.setGhost(true);
    return atom;
}

void AtomManager::removeGhostAtoms()
{
    m_ghostAtoms.removeAllAtoms();
}

void AtomManager::removeAllAtoms()
{
    m_atoms.removeAllAtoms();
    m_ghostAtoms.removeAllAtoms();
}

void AtomManager::setCutoffDistance(double cutoffDistance) {
    if(m_cellData.cutoffDistance != cutoffDistance) {
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

AtomList &AtomManager::ghostAtoms()
{
    return m_ghostAtoms;
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

    int numCellsTotal = m_cellData.numberOfCellsWithGhostCells[0]*m_cellData.numberOfCellsWithGhostCells[1]*m_cellData.numberOfCellsWithGhostCells[2];

    if(m_cellData.cells.size() < numCellsTotal) {
        m_cellData.cells.resize(numCellsTotal);
    }

    for(int cellIndex=0; cellIndex<m_cellData.cells.size(); cellIndex++) {
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
    m_cellDataDirty = false;
    if(m_cellStructureDirty) updateCellStructure();

    for(Cell &cell : m_cellData.cells) {
        cell.reset();
    }

    CellData &cellData = m_cellData;

    m_atoms.iterate([&](Atom &atom, const int &atomIndex) {
        int cellIndex = Cell::cellIndexForAtom(atom, cellData);
        Cell &cell = cellData.cells.at(cellIndex);
        cell.addAtom(&atom);
    });

    m_ghostAtoms.iterate([&](Atom &atom, const int &atomIndex) {
        int cellIndex = Cell::cellIndexForAtom(atom, cellData);
        Cell &cell = cellData.cells.at(cellIndex);
        cell.addAtom(&atom);
    });

}
