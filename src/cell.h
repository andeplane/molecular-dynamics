#ifndef SYSTEMCELL_H
#define SYSTEMCELL_H
#include <vector>
#include <atom.h>
#include <iostream>
using std::vector; using std::cout; using std::endl;
class Cell;

struct CellData {
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
public:
    Cell();
    void addAtom(Atom *atom);
    vector<Atom *> atoms() const;
    void setAtoms(const vector<Atom *> &atoms);
    void reset();

    static int cellIndexFromIJK(const int &i, const int &j, const int &k, CellData &cellData) {
        return i*cellData.numberOfCellsWithGhostCells[1]*cellData.numberOfCellsWithGhostCells[2]+j*cellData.numberOfCellsWithGhostCells[2]+k;
    }

    static int cellIndexForAtom(Atom *atom, CellData &cellData) {
        int i = atom->position[0] / cellData.cellLength[0] + 1;
        int j = atom->position[1] / cellData.cellLength[1] + 1;
        int k = atom->position[2] / cellData.cellLength[2] + 1;
        return Cell::cellIndexFromIJK(i,j,k,cellData);
    }
};

#endif // SYSTEMCELL_H
