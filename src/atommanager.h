#ifndef ATOMLIST_H
#define ATOMLIST_H

#include <vector>
#include <map>
using std::vector; using std::map; using std::pair;

#include <cell.h>
class Atom;

/*
 static int cellIndexFromIJK(const int &i, const int &j, const int &k) {
        return i*numberOfCellsWithGhostCells[1]*numberOfCellsWithGhostCells[2]+j*numberOfCellsWithGhostCells[2]+k;
    }

    static int cellIndexForAtom(Atom *atom) {
        return Cell::cellIndexForAtom(*atom);
    }

    static int cellIndexForAtom(Atom &atom) {
        int i = atom.position[0] / cellLength[0] + 1;
        int j = atom.position[1] / cellLength[1] + 1;
        int k = atom.position[2] / cellLength[2] + 1;
        cout << &atom << endl;
        cout << " has cell index: " << i << " " << j << " " << k << endl;

        return cellIndexFromIJK(i,j,k);
    }
    */

class AtomManager
{
private:
    vector<Atom> m_allAtoms;    // All atoms
    vector<Atom*> m_atoms;      // Atoms in use
    CellData m_cellData;
    map<Atom*, pair<int,int> > m_indexMap;
    bool m_cellStructureDirty;
    bool m_cellDataDirty;

    int m_nextFreeIndex;
    int indexInAtoms(Atom *atom);
    int indexInAllAtoms(Atom *atom);
    void increaseNumberOfAtoms();
    void updateCellStructure();
public:
    AtomManager(int initialAtomCount = 1000);

    Atom *addAtom();
    void removeAtom(Atom *atom);
    int atomCapacity();
    void removeAllAtoms();
    vector<Atom *> &atoms();
    void setCutoffDistance(const double &cutoffDistance);
    void setSystemLength(const vector<double> &systemLength);
    void updateCellList();
    void moveAtoms(const double &timestep);
    CellData cellData();
};

#endif // ATOMLIST_H
