#include <atommanager.h>
#include <atom.h>
#include <iostream>
#include <cmath>

using std::cerr; using std::endl; using std::cout;

AtomManager::AtomManager(int initialAtomCount) :
    m_nextFreeIndex(0),
    m_cellStructureDirty(true),
    m_cellDataDirty(true)
{
    m_allAtoms.resize(initialAtomCount);
    m_cellData.cutoffDistance = INFINITY;
    m_cellData.initialized = false;
}

AtomManager::~AtomManager()
{
    m_allAtoms.clear();
    m_atoms.clear();
    m_indexMap.clear();
}

vector<Atom *> &AtomManager::atoms()
{
    return m_atoms;
}

int AtomManager::indexInAllAtoms(Atom *atom)
{
    pair<int,int> &indexPair = m_indexMap.find(atom)->second;
    return indexPair.first;
}

void AtomManager::increaseNumberOfAtoms()
{
    int newSize = 2*(m_allAtoms.size()+1);

    vector<int> indicesOfAtomsInUse;
    indicesOfAtomsInUse.reserve(m_atoms.size());
    for(Atom *atom : m_atoms) {
        indicesOfAtomsInUse.push_back(indexInAllAtoms(atom));
    }

    m_allAtoms.resize(newSize);
    m_atoms.clear();
    m_indexMap.clear();

    for(int indexInAllAtoms : indicesOfAtomsInUse) {
        m_indexMap.insert(pair<Atom*,pair<int,int>>(&m_allAtoms.at(indexInAllAtoms), pair<int,int>(indexInAllAtoms, m_atoms.size())));
        m_atoms.push_back(&m_allAtoms.at(indexInAllAtoms));
    }

    indicesOfAtomsInUse.clear();
}


CellData &AtomManager::cellData()
{
    if(m_cellDataDirty) updateCellList();
    return m_cellData;
}

int AtomManager::indexInAtoms(Atom *atom)
{
    pair<int,int> &indexPair = m_indexMap.find(atom)->second;
    return indexPair.second;
}

Atom *AtomManager::addAtom()
{
    if(m_nextFreeIndex>=m_allAtoms.size()) {
        increaseNumberOfAtoms();
    }

    int indexInAtoms = m_atoms.size();
    int indexInAllAtoms = m_nextFreeIndex;

    Atom *nextAtom = &m_allAtoms.at(indexInAllAtoms);
    nextAtom->setMoved(false);
    nextAtom->setGhost(false);

    m_indexMap.insert(pair<Atom*,pair<int,int>>(nextAtom, pair<int,int>(indexInAllAtoms, indexInAtoms)));
    m_atoms.push_back(nextAtom);
    m_nextFreeIndex++;

    return nextAtom;
}

void AtomManager::removeAtom(Atom *atom)
{
    if(m_indexMap.count(atom) == 0) {
        cerr << "Tried to remove atom not in used. Check your code." << endl;
        return;
    }

    int removedAtomIndexAtoms = indexInAtoms(atom);
    int removedAtomIndexAllAtoms = indexInAllAtoms(atom);
    int lastUsedAtomIndexAllAtoms = --m_nextFreeIndex;

    std::swap(m_allAtoms.at(removedAtomIndexAllAtoms), m_allAtoms.at(lastUsedAtomIndexAllAtoms));
    m_atoms.pop_back();

    // allAtoms.at(removedAtomIndexAllAtoms) = allAtoms.at(lastUsedAtomIndexAllAtoms);
    Atom *removedAtomAllAtoms = &m_allAtoms.at(removedAtomIndexAllAtoms);
    m_indexMap.insert(pair<Atom*,pair<int,int>>(removedAtomAllAtoms, pair<int,int>(removedAtomIndexAllAtoms, removedAtomIndexAtoms)));

    Atom *lastUsedAtomAllAtoms = &m_allAtoms.at(lastUsedAtomIndexAllAtoms);
    m_indexMap.erase(m_indexMap.find(lastUsedAtomAllAtoms)); // Remove the pointer of the atom that is not used anymore
}

int AtomManager::atomCapacity()
{
    return m_allAtoms.size();
}

void AtomManager::removeAllAtoms()
{
    m_indexMap.clear();
    m_atoms.clear();
    m_nextFreeIndex = 0;
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

    for(Atom *atom : m_atoms) {
        int cellIndex = Cell::cellIndexForAtom(atom, m_cellData);
        Cell &cell = m_cellData.cells.at(cellIndex);
        cell.addAtom(atom);
    }
}

void AtomManager::moveAtoms(const double &timestep) {
    for(Atom *atom : m_atoms) {
        atom->move(timestep);
    }

    m_cellDataDirty = true;
}
