#include "topology.h"

#include <system.h>
#include <atom.h>
#include <algorithm>
#include <utility>
#include <atomtype.h>
#include <cmath>
#include <atommanager.h>
#include <includes.h>

Topology::Topology() :
    m_isInitialized(false),
    m_mpiReceiveBuffer(0),
    m_mpiSendBuffer(0)
{

}

Topology::~Topology()
{
    m_nodeIndices.clear();
    m_nodeLength.clear();
    m_systemLength.clear();
    for(vector<unsigned long> vec : m_moveQueue) {
         vec.clear();
    }
    m_moveQueue.clear();
    m_mpiSendBuffer.clear();
    m_mpiReceiveBuffer.clear();
}

void Topology::initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength)
{
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/
    m_mpiReceiveBuffer.resize(1e6,0);
    m_mpiSendBuffer.resize(1e6,0);
    m_nodeIndex = nodeIndex;
    m_nodeLength.resize(3);
    m_systemLength.resize(3);
    m_nodeIndices.resize(3);
    m_origo.resize(3);
    m_isInitialized = true;

    for(int a=0; a<3; a++) {
        m_numNodesVector[a] = numNodesVector[a];
        m_systemLength[a] = systemLength[a];
        m_nodeLength[a] = systemLength[a] / numNodesVector[a];
    }

    m_numNodes = numNodesVector[0]*numNodesVector[1]*numNodesVector[2];
    m_nodeIndices[0] = nodeIndex/(numNodesVector[1]*numNodesVector[2]);
    m_nodeIndices[1] = (nodeIndex/numNodesVector[2]) % numNodesVector[1];
    m_nodeIndices[2] = nodeIndex%numNodesVector[2];

    m_origo[0] = m_nodeIndices[0]*m_nodeLength[0];
    m_origo[1] = m_nodeIndices[1]*m_nodeLength[1];
    m_origo[2] = m_nodeIndices[2]*m_nodeLength[2];

    /* Integer vectors to specify the six neighbor nodes */
    int integerVector[6][3] = {
        {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
    };

    int k1[3];

    for (int n=0; n<6; n++) {
        /* Vector index of neighbor */
        for (int a=0; a<3; a++) {
            k1[a] = (m_nodeIndices[a]+integerVector[n][a]+m_numNodesVector[a])%m_numNodesVector[a];
        }

        /* Scalar neighbor ID */
        m_neighborNodesIndices[n] = k1[0]*m_numNodesVector[1]*m_numNodesVector[2]+k1[1]*m_numNodesVector[2]+k1[2];

        /* Shift vector */
        for (int a=0; a<3; a++) m_shiftVector[n][a] = m_nodeLength[a]*integerVector[n][a];
    }

    /* Set up the node parity table, myparity */
    for (int a=0; a<3; a++) {
        if (m_numNodesVector[a] == 1) {
            m_parity[a] = 2;
        } else if (m_nodeIndices[a]%2 == 0) {
            m_parity[a] = 0;
        } else {
            m_parity[a] = 1;
        }
    }

    m_moveQueue.resize(6);
}


int Topology::numNodes() const
{
    return m_numNodes;
}

int Topology::nodeIndex() const
{
    return m_nodeIndex;
}

vector<double> Topology::systemLength() const
{
    return m_systemLength;
}

void Topology::setSystemLength(vector<double> systemLength)
{
    vector<int> numNodesVector(3,0);
    numNodesVector[0] = m_numNodesVector[0]; numNodesVector[1] = m_numNodesVector[1]; numNodesVector[2] = m_numNodesVector[2];
    initialize(m_nodeIndex, numNodesVector,systemLength);
    numNodesVector.clear();
}

vector<double> Topology::nodeLength() const
{
    return m_nodeLength;
}

vector<int> Topology::nodeIndices() const
{
    return m_nodeIndices;
}


double Topology::maxSystemLength() const
{
    return std::max(m_systemLength[0],std::max(m_systemLength[1], m_systemLength[2]));
}


vector<double> Topology::origo() const
{
    return m_origo;
}

bool Topology::atomShouldBeCopied(Atom &atom, int &dimension, bool higher, double &cutoffDistance) {
    if (higher) return atom.position[dimension] > m_nodeLength[dimension]-cutoffDistance;
    else return atom.position[dimension] < cutoffDistance;
}

bool Topology::atomDidChangeNode(Atom &atom, int &dimension, bool higher) {
    if (higher) return atom.position[dimension] >= m_nodeLength[dimension];
    else return atom.position[dimension] < 0.0;
}

void Topology::MPIMove(System &system) {
    for(vector<unsigned long> &queue : m_moveQueue) {
        queue.clear();
    }

    for(int dimension=0;dimension<3;dimension++) {
        system.atomManager().atoms().iterate([&](Atom &atom, const int &atomIndex) {
            int nodeLower = 2*dimension;
            int nodeHigher = 2*dimension+1;

            if (!atom.removed()) { /* Don't scan moved-out atoms */
                if(atomDidChangeNode(atom,dimension,false)) {
                    m_moveQueue.at(nodeLower).push_back(atom.uniqueId());
                } else if (atomDidChangeNode(atom, dimension, true)) {
                    m_moveQueue.at(nodeHigher).push_back(atom.uniqueId());
                }
            }
        });

        for (int higher=0; higher<=1; higher++) {
            int localNodeID=2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();
            int numberToReceive = numberToSend;

            int i = 0;
            for (unsigned long uniqueId : m_moveQueue.at(localNodeID)) {
                Atom *atom = system.atomManager().getAtomByUniqueId(uniqueId);

                /* Shift the coordinate origin */
                m_mpiSendBuffer[12*i    + 0] = atom->position[0] - m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[12*i    + 1] = atom->position[1] - m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[12*i    + 2] = atom->position[2] - m_shiftVector[localNodeID][2];
                m_mpiSendBuffer[12*i+ 3 + 0] = atom->velocity[0];
                m_mpiSendBuffer[12*i+ 3 + 1] = atom->velocity[1];
                m_mpiSendBuffer[12*i+ 3 + 2] = atom->velocity[2];
                m_mpiSendBuffer[12*i+ 6 + 0] = atom->initial_position[0];
                m_mpiSendBuffer[12*i+ 6 + 1] = atom->initial_position[1];
                m_mpiSendBuffer[12*i+ 6 + 2] = atom->initial_position[2];
                m_mpiSendBuffer[12*i + 9]    = (double)atom->id();
                m_mpiSendBuffer[12*i + 10]   = (double)atom->type()->atomicNumber();
                m_mpiSendBuffer[12*i + 11]   = (double)atom->originalUniqueId();
                atom->setRemoved(true);
                i++;
            }

            memcpy(&m_mpiReceiveBuffer.front(),&m_mpiSendBuffer.front(),12*numberToReceive*sizeof(double));

            for (i=0; i<numberToReceive; i++) {
                Atom &atom = system.addAtom();
                atom.position[0] = m_mpiReceiveBuffer[12*i + 0];
                atom.position[1] = m_mpiReceiveBuffer[12*i + 1];
                atom.position[2] = m_mpiReceiveBuffer[12*i + 2];

                atom.velocity[0] = m_mpiReceiveBuffer[12*i + 3];
                atom.velocity[1] = m_mpiReceiveBuffer[12*i + 4];
                atom.velocity[2] = m_mpiReceiveBuffer[12*i + 5];

                atom.initial_position[0] = m_mpiReceiveBuffer[12*i + 6];
                atom.initial_position[1] = m_mpiReceiveBuffer[12*i + 7];
                atom.initial_position[2] = m_mpiReceiveBuffer[12*i + 8];

                atom.setId(m_mpiReceiveBuffer[12*i + 9]);
                int atomicNumber = m_mpiReceiveBuffer[12*i + 10];
                unsigned long originalUniqueId = m_mpiReceiveBuffer[12*i + 11];
                atom.setType(AtomType::atomTypeFromAtomicNumber(atomicNumber));
                atom.setOriginalUniqueId(originalUniqueId);
            }
        }
    }

    system.atomManager().atoms().cleanupList(); // Remove holes (removed atoms) by moving the last atoms
}

void Topology::copyGhostAtomsWithMPI(AtomManager &atomManager)
{
    atomManager.removeGhostAtoms();
    double cutoffDistance = atomManager.cellData().cutoffDistance;

    for(int dimension=0;dimension<3;dimension++) {
        for(vector<unsigned long> &queue : m_moveQueue) {
            queue.clear();
        }

        atomManager.ghostAtoms().iterate([&](Atom &atom) {
            for(int higher=0; higher<=1; higher++) {
                int localNodeID = 2*dimension + higher;
                if(atomShouldBeCopied(atom,dimension,higher,cutoffDistance)) {
                    m_moveQueue.at(localNodeID).push_back(atom.uniqueId());
                }
            }
        });

        atomManager.atoms().iterate([&](Atom &atom) {
            for(int higher=0; higher<=1; higher++) {
                int localNodeID = 2*dimension + higher;
                if(atomShouldBeCopied(atom,dimension,higher,cutoffDistance)) {
                    m_moveQueue.at(localNodeID).push_back(atom.uniqueId());
                }
            }
        });

        /* Loop through higher and lower node in this dimension */
        for(int higher=0;higher<=1;higher++) {
            int localNodeID = 2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();
            int numberToRecieve = numberToSend;

            int i = 0;
            int numberOfValuesToCopy = 6;

            for (unsigned long uniqueId : m_moveQueue.at(localNodeID)) {
                Atom *atom = atomManager.getAtomByUniqueId(uniqueId);

                m_mpiSendBuffer[numberOfValuesToCopy*i+0] = atom->position[0]-m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[numberOfValuesToCopy*i+1] = atom->position[1]-m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[numberOfValuesToCopy*i+2] = atom->position[2]-m_shiftVector[localNodeID][2];
                m_mpiSendBuffer[numberOfValuesToCopy*i+3] = atom->type()->atomicNumber();
                m_mpiSendBuffer[numberOfValuesToCopy*i+4] = atom->id();
                m_mpiSendBuffer[numberOfValuesToCopy*i+5] = atom->originalUniqueId();
                i++;
            }
            memcpy(&m_mpiReceiveBuffer.front(),&m_mpiSendBuffer.front(),numberOfValuesToCopy*numberToRecieve*sizeof(double));

            for(int i=0; i<numberToRecieve; i++) {
                Atom &ghostAtom = atomManager.addGhostAtom();
                ghostAtom.position[0] = m_mpiReceiveBuffer[numberOfValuesToCopy*i+0];
                ghostAtom.position[1] = m_mpiReceiveBuffer[numberOfValuesToCopy*i+1];
                ghostAtom.position[2] = m_mpiReceiveBuffer[numberOfValuesToCopy*i+2];
                int atomicNumber = m_mpiReceiveBuffer[numberOfValuesToCopy*i+3];
                int atomId = m_mpiReceiveBuffer[numberOfValuesToCopy*i+4];
                unsigned long OriginalUniqueId = m_mpiReceiveBuffer[numberOfValuesToCopy*i+5];

                ghostAtom.setType(AtomType::atomTypeFromAtomicNumber(atomicNumber));
                ghostAtom.setId(atomId);
                ghostAtom.setOriginalUniqueId(OriginalUniqueId);
            }
        }
    }
}
