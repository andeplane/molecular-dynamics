#ifndef ATOMMANAGER_H
#define ATOMMANAGER_H

#include <vector>
#include <map>
using std::vector; using std::map; using std::pair;

#include <atomlist.h>
#include <cell.h>
class Atom;

class AtomManager
{
private:
    AtomList m_atoms;
    AtomList m_ghostAtoms;
    CellData m_cellData;
    bool m_cellStructureDirty;
    bool m_cellDataDirty;

    void updateCellStructure();
public:
    AtomManager();
    ~AtomManager();
    Atom &addAtom();
    Atom &addGhostAtom();
    AtomList &atoms();
    AtomList &ghostAtoms();
    void removeAllAtoms();
    void removeGhostAtoms();
    void updateCellList();

    CellData &cellData();
    void setCutoffDistance(double cutoffDistance);
    double cutoffDistance();
    void setSystemLength(vector<double> &systemLength);
    int numberOfAtoms();
    int numberOfGhostAtoms();
};

#endif // ATOMMANAGER_H
