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
    static int numberOfCells[3];
    static double cellLength[3];
    static int cellIndexForAtom(Atom &atom) {
        int i = atom.position[0] / cellLength[0];
        int j = atom.position[1] / cellLength[1];
        int k = atom.position[2] / cellLength[2];

        return i*numberOfCells[1]*numberOfCells[2]+j*numberOfCells[2]+k;
    }

    SystemCell();
    void addAtom(Atom *atom);
    vector<Atom *> atoms() const;
    void setAtoms(const vector<Atom *> &atoms);
    void reset();
};

#endif // SYSTEMCELL_H
