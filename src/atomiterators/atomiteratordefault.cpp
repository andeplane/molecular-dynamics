#include "atomiteratordefault.h"
#include <cell.h>
#include <atommanager.h>
#include <includes.h>

AtomIteratorDefault::AtomIteratorDefault()
{
}

void AtomIteratorDefault::iterate(AtomManager &atomManager) {
    CellData &cellData = atomManager.cellData();
    vector<Cell> &cells = cellData.cells;

    for(int cellX=1; cellX<=cellData.numberOfCellsWithoutGhostCells[0]; cellX++) {
        for(int cellY=1; cellY<=cellData.numberOfCellsWithoutGhostCells[1]; cellY++) {
            for(int cellZ=1; cellZ<=cellData.numberOfCellsWithoutGhostCells[2]; cellZ++) {
                Cell &cell1 = cells.at(Cell::cellIndexFromIJK(cellX, cellY, cellZ, cellData));
                vector<Atom*> &atomsInCell1 = cell1.atoms();

                for(int cell2X=cellX-1; cell2X<=cellX+1; cell2X++) {
                    for(int cell2Y=cellY-1; cell2Y<=cellY+1; cell2Y++) {
                        for(int cell2Z=cellZ-1; cell2Z<=cellZ+1; cell2Z++) {
                            Cell &cell2 = cells.at(Cell::cellIndexFromIJK(cell2X, cell2Y, cell2Z, cellData));
                            vector<Atom*> &atomsInCell2 = cell2.atoms();

                            for(unsigned long atomIndex1 = 0; atomIndex1<atomsInCell1.size(); atomIndex1++) {
                                Atom *atom1 = safeOrQuickVectorLookup(atomsInCell1, atomIndex1);
                                for(unsigned long atomIndex2 = 0; atomIndex2<atomsInCell2.size(); atomIndex2++) {
                                    Atom *atom2 = safeOrQuickVectorLookup(atomsInCell2, atomIndex2);
                                    if(atom1->originalUniqueId() <= atom2->originalUniqueId() && !atom2->ghost()) continue; // Newton's 3rd law, always calculate if atom2 is ghost
                                    this->twoParticleAction()(atom1,atom2);
                                }
                            } // Loop atoms
                        }
                    }
                } // Loop neighbor cells
            }
        }
    }

    // TODO: Add 3-particle loops as well
}
