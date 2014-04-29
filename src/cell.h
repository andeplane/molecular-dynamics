#ifndef SYSTEMCELL_H
#define SYSTEMCELL_H
#include <vector>
#include <atom.h>
#include <iostream>
using std::vector; using std::cout; using std::endl;
class Cell;

struct CellData {
    bool initialized;
    vector<Cell> cells;
    int numberOfCellsWithoutGhostCells[3];
    int numberOfCellsWithGhostCells[3];
    double cellLength[3];
    double systemLength[3];
    double cutoffDistance;
};

class Cell
{
private:
    vector<Atom*> m_atoms;
    int m_cellIndex;
public:
    Cell();
    void addAtom(Atom *atom);
    vector<Atom *> atoms() const;
    void setAtoms(const vector<Atom *> &atoms);
    void reset();

    static int cellIndexFromIJK(const int &i, const int &j, const int &k, CellData &cellData) {
        int index = i*cellData.numberOfCellsWithGhostCells[1]*cellData.numberOfCellsWithGhostCells[2]+j*cellData.numberOfCellsWithGhostCells[2]+k;
        return index;
    }

    static int cellIndexForAtom(Atom *atom, CellData &cellData) {
        int i = atom->position[0] / cellData.cellLength[0] + 1;
        int j = atom->position[1] / cellData.cellLength[1] + 1;
        int k = atom->position[2] / cellData.cellLength[2] + 1;
        return Cell::cellIndexFromIJK(i,j,k,cellData);
    }

    static Cell *cellContainingAtom(Atom *atom, CellData &cellData) {
        for(Cell &cell : cellData.cells) {
            if(std::find(cell.atoms().begin(),cell.atoms().end(), atom) != cell.atoms().end()) {
                return &cell;
            }
        }

        cout << "There are currently no cells containing " << *atom << endl;
        return 0;
    }

    int cellIndex() const;
    void setCellIndex(int cellIndex);
};

#endif // SYSTEMCELL_H
