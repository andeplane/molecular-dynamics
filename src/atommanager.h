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
    vector<Atom> m_atoms;    // All atoms
    vector<Atom> m_ghostAtoms;    // All atoms
    CellData m_cellData;
    map<Atom*, pair<int,int> > m_indexMap;
    bool m_cellStructureDirty;
    bool m_cellDataDirty;

    int indexInAtoms(Atom *atom);
    int indexInAllAtoms(Atom *atom);
    void increaseNumberOfAtoms();
    void updateCellStructure();
    Atom *getNextAtom();
public:
    int m_nextFreeIndex;
    AtomManager(int initialAtomCount = 1000);
    ~AtomManager();
    vector<Atom *> &atoms();
    vector<Atom *> &ghostAtoms();
    Atom *addAtom();
    Atom *addGhostAtom();
    void removeAtom(Atom *atom);
    void removeGhostAtoms();

    int atomCapacity();
    void removeAllAtoms();
    void setCutoffDistance(double cutoffDistance);
    double cutoffDistance();
    void setSystemLength(vector<double> &systemLength);
    void updateCellList();
    void moveAtoms(const double &timestep);
    CellData &cellData();
};

#endif // ATOMLIST_H
