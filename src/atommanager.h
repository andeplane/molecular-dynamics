#pragma once
#include <map>
using std::map; using std::pair;

#include <atomlist.h>
#include <cell.h>
class Atom; class Topology;

class AtomManager
{
private:
    CellData m_cellData;
    AtomList m_atoms;
    AtomList m_ghostAtoms;
    friend std::ostream& operator<<(std::ostream&stream, AtomManager &atomManager);
    bool m_cellStructureDirty;
    bool m_cellDataDirty;
    bool m_ghostAtomsDirty;
    bool m_updatingGhostAtoms;
    bool m_ghostAtomsEnabled;
    function<void(Atom &atom)> m_onAtomAdd;
    vector<double> m_mpiSendBuffer;
    vector<double> m_mpiReceiveBuffer;
    Topology *m_topology;


    void updateCellStructure();
public:
    AtomManager();
    ~AtomManager();
    Atom &addAtom(shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    Atom &addGhostAtom(shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::NoAtom));
    AtomList &atoms();
    AtomList &ghostAtoms();
    void removeAllAtoms();
    void removeGhostAtoms();
    void updateCellList();
    void applyPeriodicBoundaryConditions(vector<double> systemLength);

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
    Atom &getAtomByUniqueId(unsigned long uniqueId);
};
