#pragma once
#include <particles/atom.h>
#include <iostream>
#include <utils/vec3.h>
using CompPhys::vec3;
using std::cout; using std::endl;

class Cell;
struct CellData {
    bool initialized;
    vector<Cell> cells;
    int numberOfCellsWithoutGhostCells[3];
    int numberOfCellsWithGhostCells[3];
    vec3 cellLength;
    vec3 systemLength;
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
    vector<Atom *> &atoms();
    void setAtoms(const vector<Atom *> &atoms);
    void reset();

    static int cellIndexFromIJK(const int &i, const int &j, const int &k, CellData &cellData) {
        int index = i*cellData.numberOfCellsWithGhostCells[1]*cellData.numberOfCellsWithGhostCells[2]+j*cellData.numberOfCellsWithGhostCells[2]+k;
        return index;
    }

    static int cellIndexForAtom(Atom &atom, CellData &cellData) {
        int i = atom.position.x() / cellData.cellLength.x() + 1;
        int j = atom.position.y() / cellData.cellLength.y() + 1;
        int k = atom.position.z() / cellData.cellLength.z() + 1;
        int index = Cell::cellIndexFromIJK(i,j,k,cellData);
        return index;
    }

    static Cell *cellContainingAtom(Atom &atom, CellData &cellData) {
        for(Cell &cell : cellData.cells) {
            if(std::find(cell.atoms().begin(),cell.atoms().end(), &atom) != cell.atoms().end()) {
                return &cell;
            }
        }

        cout << "There are currently no cells containing " << atom << endl;
        return 0;
    }

    int cellIndex() const;
    void setCellIndex(int cellIndex);
};
