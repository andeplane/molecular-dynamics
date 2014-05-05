#ifndef ATOMMANAGER_H
#define ATOMMANAGER_H

#include <vector>
#include <map>
using std::vector; using std::map; using std::pair;

#include <atomlist.h>
#include <cell.h>
class Atom; class Topology;

class AtomManager
{
private:
    friend std::ostream& operator<<(std::ostream&stream, AtomManager &atomManager);
    AtomList m_atoms;
    AtomList m_ghostAtoms;
    CellData m_cellData;
    bool m_cellStructureDirty;
    bool m_cellDataDirty;
    bool m_ghostAtomsDirty;
    bool m_updatingGhostAtoms;
    bool m_ghostAtomsEnabled;
    vector<double> m_mpiSendBuffer;
    vector<double> m_mpiReceiveBuffer;
    Topology *m_topology;

    void updateCellStructure();
public:
    AtomManager();
    ~AtomManager();
    Atom &addAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    Atom &addGhostAtom(AtomType *atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
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
    void setTopology(Topology *topology);
    void updateGhostAtoms();
    bool ghostAtomsEnabled() const;
    void setGhostAtomsEnabled(bool ghostAtomsEnabled);
    Atom *getAtomByOriginalUniqueId(unsigned long uniqueId);
    Atom *getOriginalAtomFromGhostAtom(Atom &ghostAtom);
    Atom *getAtomByUniqueId(unsigned long uniqueId);
};

#endif // ATOMMANAGER_H
