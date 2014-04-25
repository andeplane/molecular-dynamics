#ifndef SYSTEMCELL_H
#define SYSTEMCELL_H
#include <vector>
#include <atom.h>

using std::vector;


class SystemCell
{
private:
    vector<Atom*> m_atoms;

public:
    static int numberOfCellsWithGhostCells[3];
    static int numberOfCellsWithoutGhostCells[3];
    static double cellLength[3];

    static int cellIndexFromIJK(const int &i, const int &j, const int &k) {
        return i*numberOfCellsWithGhostCells[1]*numberOfCellsWithGhostCells[2]+j*numberOfCellsWithGhostCells[2]+k;
    }

    static int cellIndexForAtom(Atom &atom) {
        int i = atom.position[0] / cellLength[0] + 1;
        int j = atom.position[1] / cellLength[1] + 1;
        int k = atom.position[2] / cellLength[2] + 1;

        return cellIndexFromIJK(i,j,k);
    }

    SystemCell();
    void addAtom(Atom *atom);
    vector<Atom *> atoms() const;
    void setAtoms(const vector<Atom *> &atoms);
    void reset();
};

#endif // SYSTEMCELL_H
