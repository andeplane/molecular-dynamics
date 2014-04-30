#include "system.h"
#include "atom.h"
#include "topology.h"
#include "cmath"
#include "potentials/lennardjonespotential.h"
#include <iostream>
using namespace std;


AtomManager &System::atomManager()
{
    return m_atomManager;
}

System::System() :
    m_indexOfFirstGhostAtom(-1),
    m_isInitialized(false)
{

}

System::~System()
{
    m_potentials.clear();

}

void System::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength)
{
    m_isInitialized = true;
    m_indexOfFirstGhostAtom = -1;
    m_topology.initialize(nodeIndex, numNodesVector, systemLength);
    m_atomManager.setSystemLength(systemLength);
}

Atom *System::addAtom()
{
    return m_atomManager.addAtom();
}

void System::removeAtom(Atom *atom)
{
    m_atomManager.removeAtom(atom);
}

void System::removeAllAtoms()
{
    m_atomManager.removeAllAtoms();
}

int System::numberOfAtoms()
{
    return m_atomManager.atoms().size();
}

int System::numberOfGhostAtoms()
{
    return m_indexOfFirstGhostAtom-numberOfAtoms();
}

void System::setSystemLength(vector<double> &systemLength)
{
    m_topology.setSystemLength(systemLength);
    m_atomManager.setSystemLength(systemLength);
}

void System::addPotential(PotentialType type) {
    checkIfInitialized();
    if(type == PotentialType::LennardJones) {
        LennardJonesPotential *lennardJones = new LennardJonesPotential();
        potentials().push_back(lennardJones);
    }
}

void System::checkIfInitialized() {
    if(!m_isInitialized) std::cerr << "Error, System object not initialized." << std::endl;
}

int System::firstGhostAtomIndex() const
{
    return m_indexOfFirstGhostAtom;
}

void System::setFirstGhostAtomIndex(int firstGhostAtomIndex)
{
    m_indexOfFirstGhostAtom = firstGhostAtomIndex;
}

vector<Potential*> &System::potentials()
{
    return m_potentials;
}

void System::resetForces() {
    checkIfInitialized();
    for(Atom *atom : m_atomManager.atoms()) {
        atom->resetForce();
    }
}

vector<Atom*> &System::atoms()
{
    return m_atomManager.atoms();
}

Topology &System::topology()
{
    return m_topology;
}
