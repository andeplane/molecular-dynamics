#include <system.h>
#include <atom.h>
#include <topology.h>
#include <cmath>
#include <potentials/lennardjonespotential.h>
#include <iostream>
using namespace std;

AtomManager &System::atomManager()
{
    return m_atomManager;
}

System::System() :
    m_isInitialized(false)
{

}

System::~System()
{
    m_potentials.clear();
    m_atomManager.removeAllAtoms();
}

void System::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength)
{
    m_isInitialized = true;
    m_topology.initialize(nodeIndex, numNodesVector, systemLength);
    m_atomManager.setTopology(&m_topology);
    m_atomManager.setSystemLength(systemLength);
    m_atomManager.removeAllAtoms();
}

Atom &System::addAtom(AtomType *atomType)
{
    checkIfInitialized();
    return m_atomManager.addAtom(atomType);
}

Atom &System::addGhostAtom(AtomType *atomType)
{
    checkIfInitialized();
    return m_atomManager.addGhostAtom(atomType);
}

void System::removeAllAtoms() {
    m_atomManager.removeAllAtoms();
}

void System::removeGhostAtoms() {
    m_atomManager.removeGhostAtoms();
}

int System::numberOfGhostAtoms()
{
    return m_atomManager.numberOfGhostAtoms();
}

int System::numberOfAtoms()
{
    return m_atomManager.numberOfAtoms();
}

void System::setSystemLength(vector<double> &systemLength)
{
    m_topology.setSystemLength(systemLength);
    m_atomManager.setSystemLength(systemLength);
}

Potential *System::addPotential(PotentialType type) {
    checkIfInitialized();
    if(type == PotentialType::LennardJones) {
        LennardJonesPotential *lennardJones = new LennardJonesPotential();
        potentials().push_back(lennardJones);
        return lennardJones;
    }
}

void System::checkIfInitialized() {
    if(!m_isInitialized) std::cerr << "Error, System object not initialized." << std::endl;
}

vector<Potential*> &System::potentials()
{
    return m_potentials;
}

Topology &System::topology()
{
    return m_topology;
}
