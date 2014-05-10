#include "atomiteratordefault.h"
#include <cell.h>
#include <atommanager.h>
#include <includes.h>

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
    CellData &cellData = atomManager.cellData();
    vector<Cell> &cells = cellData.cells;
    double maximumNeighborDistanceSquared = m_maximumNeighborDistance*m_maximumNeighborDistance;

    cout << "Starting to iterate..." << endl;

    for(int cellX=0; cellX<cellData.numberOfCellsWithGhostCells[0]; cellX++) {
        for(int cellY=0; cellY<cellData.numberOfCellsWithGhostCells[1]; cellY++) {
            for(int cellZ=0; cellZ<cellData.numberOfCellsWithGhostCells[2]; cellZ++) {
                Cell &cell1 = safeOrQuickVectorLookup(cells,Cell::cellIndexFromIJK(cellX, cellY, cellZ, cellData));

                for(int cell2X=cellX-1; cell2X<=cellX+1; cell2X++) {
                    for(int cell2Y=cellY-1; cell2Y<=cellY+1; cell2Y++) {
                        for(int cell2Z=cellZ-1; cell2Z<=cellZ+1; cell2Z++) {

                            int cell2XPeriodic = (cell2X + cellData.numberOfCellsWithGhostCells[0]) % cellData.numberOfCellsWithGhostCells[0];
                            int cell2YPeriodic = (cell2Y + cellData.numberOfCellsWithGhostCells[1]) % cellData.numberOfCellsWithGhostCells[1];
                            int cell2ZPeriodic = (cell2Z + cellData.numberOfCellsWithGhostCells[2]) % cellData.numberOfCellsWithGhostCells[2];
                            Cell &cell2 = safeOrQuickVectorLookup(cells,Cell::cellIndexFromIJK(cell2XPeriodic, cell2YPeriodic, cell2ZPeriodic, cellData));

                            for(Atom *atom1 : cell1.atoms()) {
                                for(Atom *atom2 : cell2.atoms()) {
                                    if(m_createNeighborList && atom1->uniqueId() < atom2->uniqueId()) {
                                        double dx = atom1->position[0] - atom2->position[0];
                                        double dy = atom1->position[1] - atom2->position[1];
                                        double dz = atom1->position[2] - atom2->position[2];
                                        double dr2 = dx*dx + dy*dy + dz*dz;

                                        if(dr2 < maximumNeighborDistanceSquared) {
                                            atom1->addNeighbor(*atom2);
                                            atom2->addNeighbor(*atom1);
                                        }
                                    }

                                    // Skip two particle forces between ghosts. Also skip two particle force between non ghosts unless atom1.uniqueId<atom2.uniqueId (the other permutation will compute forces)
                                    if( (atom1->ghost() && atom2->ghost()) || (atom1->uniqueId() >= atom2->uniqueId() && (!atom1->ghost() && !atom2->ghost()) )) continue;
                                    m_twoParticleAction(atom1,atom2);
                                }
                            } // Loop atoms
                        }
                    }
                } // Loop neighbor cells
            }
        }
    }

    if(!m_threeParticleAction) return;

    // Three particle loop
    auto threeParticleLooper = [&](Atom &atom) {
        Atom *atom1 = &atom;
        vector<atomUniqueId> &neighbors = atom1->neighbors();
        if(neighbors.size() < 2) return; // No three particle contribution here

        atomUniqueId atom1UniqueId = atom1->uniqueId();
        for(unsigned long neighborIndex1 = 0; neighborIndex1<neighbors.size(); neighborIndex1++) {
            atomUniqueId &atom2UniqueId = safeOrQuickVectorLookup(neighbors, neighborIndex1);
            Atom *atom2 = &atomManager.getAtomByUniqueId(atom2UniqueId);
            if(atom1UniqueId >= atom2UniqueId) continue; // Only accept configurations where the unique atom ids are in ascending order

            for(unsigned long neighborIndex2 = 0; neighborIndex2<neighbors.size(); neighborIndex2++) {
                atomUniqueId &atom3UniqueId = safeOrQuickVectorLookup(neighbors, neighborIndex2);
                Atom *atom3 = &atomManager.getAtomByUniqueId(atom3UniqueId);

                if(atom2UniqueId >= atom3UniqueId) continue;  // Only accept configurations where the unique atom ids are in ascending order
                if(atom1->ghost() && atom2->ghost() && atom3->ghost()) continue;

                m_threeParticleAction(atom1,atom2,atom3);
            }
        }
    };

    cout << "Looping through atoms..." << endl;
    atomManager.atoms().iterate(threeParticleLooper);
    cout << "Looping through ghost atoms..." << endl;
    atomManager.ghostAtoms().iterate(threeParticleLooper);
}
