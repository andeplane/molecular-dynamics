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
//    cout << endl << endl << endl << "Will iterate through all atoms. List:" << atomManager << endl << endl;
//    cout << "Max neighbor distance: " << UnitConverter::lengthToAngstroms(m_maximumNeighborDistance) << endl;

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
                                    // cout << "Pair: (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << "," << atom2->uniqueId() << (atom2->ghost() ? "*" : "") << ")" << endl;
                                    if(m_createNeighborList && atom1->uniqueId() < atom2->uniqueId()) {
                                        double dx = atom1->position[0] - atom2->position[0];
                                        double dy = atom1->position[1] - atom2->position[1];
                                        double dz = atom1->position[2] - atom2->position[2];
                                        double dr2 = dx*dx + dy*dy + dz*dz;

                                        if(dr2 < maximumNeighborDistanceSquared) {
                                            // cout << "Added pair: (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << "," << atom2->uniqueId() << (atom2->ghost() ? "*" : "") << ")" << endl;
                                            atom1->addNeighbor(*atom2);
                                            atom2->addNeighbor(*atom1);
                                        }
                                    }

                                    // Skip two particle forces between ghosts. Also skip two particle force between non ghosts unless atom1.uniqueId<atom2.uniqueId (the other permutation will compute forces)
                                    if( (atom1->ghost() && atom2->ghost()) || (atom1->uniqueId() <= atom2->uniqueId() && (!atom1->ghost() && !atom2->ghost()) )) continue;
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

//    atomManager.atoms().iterate([] (Atom &atom) {
//       cout << "Neighbors " << atom.uniqueId() << ": " << atom.neighbors() << endl;
//    });

    auto threeParticleLooper = [&](Atom &atom) {
        Atom *atom1 = &atom;
//        cout << endl << endl << "Num neighbors = " << atom1->neighbors().size() << " for " << atom << endl;
//        cout << "Neighbors: " << atom1->neighbors() << endl;
        vector<atomUniqueId> &neighbors = atom1->neighbors();
        if(neighbors.size() < 2) return; // No three particle contribution here

        atomUniqueId atom1UniqueId = atom1->uniqueId();
        // cout << "Singlet: " << "(" << atom1OriginalUniqueId << (atom1->ghost() ? "*" : "") << ")  --- unique id (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << ")" << endl;

        for(unsigned long neighborIndex1 = 0; neighborIndex1<neighbors.size()-1; neighborIndex1++) {
            atomUniqueId &atom2UniqueId = safeOrQuickVectorLookup(neighbors, neighborIndex1);
            Atom *atom2 = &atomManager.getAtomByUniqueId(atom2UniqueId);

            // cout << "Pair: (" << atom1OriginalUniqueId << (atom1->ghost() ? "*" : "") << "," << atom2OriginalUniqueId << (atom2->ghost() ? "*" : "") << ")" << "  --- unique ids (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << "," << atom2->uniqueId() << (atom2->ghost() ? "*" : "") << ")" << endl;
            if(atom1UniqueId >= atom2UniqueId) continue; // Only accept configurations where the original atom ids are in ascending order

            for(unsigned long neighborIndex2 = neighborIndex1+1; neighborIndex2<neighbors.size(); neighborIndex2++) {
                atomUniqueId &atom3UniqueId = safeOrQuickVectorLookup(neighbors, neighborIndex2);
                Atom *atom3 = &atomManager.getAtomByUniqueId(atom3UniqueId);

                // cout << "Triplet: (" << atom1OriginalUniqueId << (atom1->ghost() ? "*" : "") << "," << atom2OriginalUniqueId << (atom2->ghost() ? "*" : "") << "," << atom3OriginalUniqueId << (atom3->ghost() ? "*" : "") << ")  --- unique ids (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << "," << atom2->uniqueId() << (atom2->ghost() ? "*" : "") << ", " << atom3->uniqueId() << (atom3->ghost() ? "*" : "") << ")" << endl;
                if(atom2UniqueId >= atom3UniqueId) continue;  // Only accept configurations where the original atom ids are in ascending order
                // cout << "Executing three particle action with ";
                // cout << " triplet: (" << atom1OriginalUniqueId << (atom1->ghost() ? "*" : "") << "," << atom2OriginalUniqueId << (atom2->ghost() ? "*" : "") << "," << atom3OriginalUniqueId << (atom3->ghost() ? "*" : "") << ")  --- unique ids (" << atom1->uniqueId() << (atom1->ghost() ? "*" : "") << "," << atom2->uniqueId() << (atom2->ghost() ? "*" : "") << ", " << atom3->uniqueId() << (atom3->ghost() ? "*" : "") << ")" << endl;
                if(atom1->ghost() && atom2->ghost() && atom3->ghost()) continue;

                m_threeParticleAction(atom1,atom2,atom3);
            }
        }
    };

    atomManager.atoms().iterate(threeParticleLooper);
    atomManager.ghostAtoms().iterate(threeParticleLooper);
    return;

    atomManager.atoms().iterate([&](Atom &atom1) {
        atomManager.atoms().iterate([&](Atom &atom2) {
            atomManager.atoms().iterate([&](Atom &atom3) {
                if(atom1.originalUniqueId() < atom2.originalUniqueId() && atom1.originalUniqueId() < atom3.uniqueId() && atom2.uniqueId() < atom3.uniqueId()) m_threeParticleAction(&atom1,&atom2,&atom3);
            });
        });

        atomManager.atoms().iterate([&](Atom &atom2) {
            atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                if(atom1.originalUniqueId() < atom2.originalUniqueId() && atom1.originalUniqueId() < atom3.uniqueId() && atom2.uniqueId() < atom3.uniqueId()) m_threeParticleAction(&atom1,&atom2,&atom3);
            });
        });

        atomManager.ghostAtoms().iterate([&](Atom &atom2) {
            atomManager.ghostAtoms().iterate([&](Atom &atom3) {
                if(atom1.originalUniqueId() < atom2.originalUniqueId() && atom1.originalUniqueId() < atom3.uniqueId() && atom2.uniqueId() < atom3.uniqueId()) m_threeParticleAction(&atom1,&atom2,&atom3);
            });
        });
    });
}
