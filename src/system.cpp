#include <system.h>
#include <atom.h>
#include <topology.h>
#include <cmath>
#include <potentials/lennardjonespotential.h>
#include <iostream>
#include <statisticssampler.h>
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
    setSystemLength(systemLength);
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

vector<double> System::systemLength()
{
    return m_topology.systemLength();
}

void System::removeTotalMomentum()
{
    StatisticsSampler sampler;
    vector<double> momentum = sampler.calculateTotalMomentum(*this);

    if(numberOfAtoms() > 0) {
        momentum[0] /= numberOfAtoms();
        momentum[1] /= numberOfAtoms();
        momentum[2] /= numberOfAtoms();
    }

    atomManager().atoms().iterate([&](Atom &atom) {
        atom.velocity[0] -= momentum[0]*atom.type()->massInverse();
        atom.velocity[1] -= momentum[1]*atom.type()->massInverse();
        atom.velocity[2] -= momentum[2]*atom.type()->massInverse();
    });
}

int System::numberOfAtoms()
{
    return m_atomManager.numberOfAtoms();
}

void System::setSystemLength(vector<double> systemLength)
{
    m_topology.setSystemLength(systemLength);
    m_atomManager.setSystemLength(systemLength);
    for(Potential *potential : potentials()) {
        potential->setSystemLength(systemLength);
    }
}

Potential *System::addPotential(PotentialType type) {
    checkIfInitialized();
    if(type == PotentialType::LennardJones) {
        LennardJonesPotential *lennardJones = new LennardJonesPotential();
        lennardJones->setSystemLength(systemLength());
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

std::ostream& operator<<(std::ostream &stream, System &system) {
    StatisticsSampler sampler;

    stream << "System information:" << endl;
    stream << "Length: (" << system.systemLength()[0] << ", " << system.systemLength()[1] << ", " << system.systemLength()[2] << ")" << endl;
    stream << "Number of atoms: " << system.numberOfAtoms() << endl;
    stream << "Number of ghost atoms: " << system.numberOfGhostAtoms() << endl;
    stream << "Potentials (" << system.potentials().size() << "): " << endl;
    for(Potential *potential : system.potentials()) {
        stream << "   " << potential->name() << endl;
    }
    vector<double> momentum = sampler.calculateTotalMomentum(system);
    stream << "Total momentum: (" << momentum.at(0) << "," << momentum.at(1) << "," << momentum.at(2) << ")";
    return stream;
}

std::ostream& operator<<(std::ostream &stream, const std::vector<double> &vec) {
    stream << "(";
    for(int i=0; i<vec.size(); i++) {
        if(i+1 == vec.size()) stream << vec[i] << ")";
        else stream << vec[i] << ", ";
    }
}
