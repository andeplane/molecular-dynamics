#include "atomiteratordefault.h"
#include <cell.h>
#include <atommanager.h>
#include <includes.h>
#include <neighborlist.h>


double AtomIteratorDefault::maximumNeighborDistance() const
{
    return m_maximumNeighborDistance;
}

void AtomIteratorDefault::setMaximumNeighborDistance(double maximumNeighborDistance)
{
    m_maximumNeighborDistance = maximumNeighborDistance;
}

bool AtomIteratorDefault::createNeighborList() const
{
    return m_createNeighborList;
}

void AtomIteratorDefault::setCreateNeighborList(bool createNeighborList)
{
    m_createNeighborList = createNeighborList;
}

void AtomIteratorDefault::setThreeParticleAction(const function<void (Atom *, Atom *, Atom *)> &threeParticleAction)
{
    AtomIterator::setThreeParticleAction(threeParticleAction);
    m_createNeighborList = true;
}

AtomIteratorDefault::AtomIteratorDefault() :
    m_maximumNeighborDistance(INFINITY),
    m_createNeighborList(false)
{

}

void AtomIteratorDefault::iterate(AtomManager &atomManager) {
    NeighborList neighborList;


    CellData &cellData = atomManager.cellData();
    vector<Cell> &cells = cellData.cells;
    double maximumNeighborDistanceSquared = m_maximumNeighborDistance*m_maximumNeighborDistance;

    for(int cellX=1; cellX<=cellData.numberOfCellsWithoutGhostCells[0]; cellX++) {
        for(int cellY=1; cellY<=cellData.numberOfCellsWithoutGhostCells[1]; cellY++) {
            for(int cellZ=1; cellZ<=cellData.numberOfCellsWithoutGhostCells[2]; cellZ++) {
                Cell &cell1 = safeOrQuickVectorLookup(cells,Cell::cellIndexFromIJK(cellX, cellY, cellZ, cellData));

                for(int cell2X=cellX-1; cell2X<=cellX+1; cell2X++) {
                    for(int cell2Y=cellY-1; cell2Y<=cellY+1; cell2Y++) {
                        for(int cell2Z=cellZ-1; cell2Z<=cellZ+1; cell2Z++) {
                            Cell &cell2 = safeOrQuickVectorLookup(cells,Cell::cellIndexFromIJK(cell2X, cell2Y, cell2Z, cellData));

                            for(Atom *atom1 : cell1.atoms()) {
                                for(Atom *atom2 : cell2.atoms()) {
                                    if(atom1->originalUniqueId() <= atom2->originalUniqueId() && !atom2->ghost()) continue; // Newton's 3rd law, always calculate if atom2 is ghost
                                    m_twoParticleAction(atom1,atom2);
                                    if(m_createNeighborList) {
                                        double dx = atom1->position[0] - atom2->position[0];
                                        double dy = atom1->position[1] - atom2->position[1];
                                        double dz = atom1->position[2] - atom2->position[2];
                                        double dr2 = dx*dx + dy*dy + dz*dz;

                                        if(dr2 < maximumNeighborDistanceSquared) {
                                            neighborList.addNeighbors(atom1,atom2);
                                        }
                                    }
                                }
                            } // Loop atoms
                        }
                    }
                } // Loop neighbor cells
            }
        }
    }
    return;
    if(!m_threeParticleAction) return;
    // Three particle loop
    atomManager.atoms().iterate([&](Atom &atom) {
        Atom *atom1 = &atom;
        vector<Atom*> neighbors = neighborList.neighborsForAtom(atom1);
        if(neighbors.size() < 2) return; // No three particle contribution here
        for(unsigned long neighborIndex1 = 0; neighborIndex1<neighbors.size()-1; neighborIndex1++) {
            Atom *atom2 = safeOrQuickVectorLookup(neighbors,neighborIndex1);
            for(unsigned long neighborIndex2 = neighborIndex1+1; neighborIndex2<neighbors.size(); neighborIndex2++) {
                Atom *atom3 = safeOrQuickVectorLookup(neighbors,neighborIndex2);
                m_threeParticleAction(atom1,atom2,atom3);
            }
        }
    });
}
