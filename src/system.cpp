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
    m_indexOfNextFreeAtom(0),
    m_firstGhostAtomIndex(-1),
    m_cutoffDistance(1),
    m_isInitialized(false)
{

}

void System::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength, double cutoffDistance)
{
    m_isInitialized = true;
    m_cutoffDistance = cutoffDistance;
    m_firstGhostAtomIndex = -1;
    m_topology.initialize(nodeIndex, numNodesVector, systemLength);
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

double System::cutoffDistance() const
{
    return m_cutoffDistance;
}

void System::setCutoffDistance(double cutoffDistance)
{
    m_cutoffDistance = cutoffDistance;
}

int System::firstGhostAtomIndex() const
{
    return m_firstGhostAtomIndex;
}

void System::setFirstGhostAtomIndex(int firstGhostAtomIndex)
{
    m_firstGhostAtomIndex = firstGhostAtomIndex;
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

Topology System::topology() const
{
    return m_topology;
}

void System::setTopology(const Topology &topology)
{
    m_topology = topology;
}
