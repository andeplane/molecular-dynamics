#include "topology.h"

#include <system.h>
#include <atom.h>
#include <algorithm>
#include <utility>
#include <atomtype.h>
#include <cmath>

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
    for(vector<Atom*> vec : m_moveQueue) {
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

bool Topology::atomShouldBeCopied(Atom &atom, int &dimension, bool higher, double &cutoffDistance) {
    if (higher) return atom.position[dimension] > m_nodeLength[dimension]-cutoffDistance;
    else return atom.position[dimension] < cutoffDistance;
}

bool Topology::atomDidChangeNode(Atom &atom, int &dimension, bool higher) {
    if (higher) return atom.position[dimension] >= m_nodeLength[dimension];
    else return atom.position[dimension] < 0.0;
}

void Topology::MPICopy(System &system, double cutoffDistance)
{
    system.removeGhostAtoms();
    cutoffDistance = std::min(cutoffDistance, maxSystemLength());

    for(vector<Atom*> &queue : m_moveQueue) {
        queue.clear();
    }

    for(int dimension=0;dimension<3;dimension++) {
        system.atomManager().atoms().iterate([&](Atom &atom, const int &atomIndex) {
            for(int higher=0; higher<=1; higher++) {
                int localNodeID = 2*dimension + higher;
                if(atomShouldBeCopied(atom,dimension,higher,cutoffDistance)) m_moveQueue.at(localNodeID).push_back(&atom);
            }
        });

        /* Loop through higher and lower node in this dimension */
        for(int higher=0;higher<2;higher++) {
            int localNodeID = 2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();
            int numberToRecieve = numberToSend;

            int i = 0;
            for (Atom *atom : m_moveQueue.at(localNodeID)) {
                /* Shift the coordinate origin */
                m_mpiSendBuffer[4*i+0] = atom->position[0]-m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[4*i+1] = atom->position[1]-m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[4*i+2] = atom->position[2]-m_shiftVector[localNodeID][2];
                m_mpiSendBuffer[4*i+3] = atom->type()->atomicNumber();
                i++;
            }

            memcpy(&m_mpiReceiveBuffer.front(),&m_mpiSendBuffer.front(),3*numberToRecieve*sizeof(double));

            for(int i=0; i<numberToRecieve; i++) {
                Atom &ghostAtom = system.addGhostAtom();
                ghostAtom.position[0] = m_mpiReceiveBuffer[4*i+0];
                ghostAtom.position[1] = m_mpiReceiveBuffer[4*i+1];
                ghostAtom.position[2] = m_mpiReceiveBuffer[4*i+2];
                int atomicNumber = m_mpiReceiveBuffer[4*i+3];
                ghostAtom.setType(AtomType::atomTypeFromAtomicNumber(atomicNumber));
            }
        }
    }
}

void Topology::MPIMove(System &system) {
    for(vector<Atom*> &queue : m_moveQueue) {
        queue.clear();
    }

    for(int dimension=0;dimension<3;dimension++) {
        system.atomManager().atoms().iterate([&](Atom &atom, const int &atomIndex) {
            int nodeLower = 2*dimension;
            int nodeHigher = 2*dimension+1;

            if (!atom.removed()) { /* Don't scan moved-out atoms */
                if(atomDidChangeNode(atom,dimension,false)) {
                    m_moveQueue.at(nodeLower).push_back(&atom);
                } else if (atomDidChangeNode(atom, dimension, true)) {
                    m_moveQueue.at(nodeHigher).push_back(&atom);
                }
            }
        });

        for (int higher=0; higher<=1; higher++) {
            int localNodeID=2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();
            int numberToReceive = numberToSend;

            /* Message buffering */
            int i = 0;
            for (Atom *atom : m_moveQueue.at(localNodeID)) {
                /* Shift the coordinate origin */
                m_mpiSendBuffer[11*i    + 0] = atom->position[0] - m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[11*i    + 1] = atom->position[1] - m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[11*i    + 2] = atom->position[2] - m_shiftVector[localNodeID][2];
                m_mpiSendBuffer[11*i+ 3 + 0] = atom->velocity[0];
                m_mpiSendBuffer[11*i+ 3 + 1] = atom->velocity[1];
                m_mpiSendBuffer[11*i+ 3 + 2] = atom->velocity[2];
                m_mpiSendBuffer[11*i+ 6 + 0] = atom->initial_position[0];
                m_mpiSendBuffer[11*i+ 6 + 1] = atom->initial_position[1];
                m_mpiSendBuffer[11*i+ 6 + 2] = atom->initial_position[2];
                m_mpiSendBuffer[11*i + 9]    = (double)atom->id();
                m_mpiSendBuffer[11*i + 10]   = (double)atom->type()->atomicNumber();
                atom->setRemoved(true);
                i++;
            }

            memcpy(&m_mpiReceiveBuffer.front(),&m_mpiSendBuffer.front(),11*numberToReceive*sizeof(double));

            /* Message storing */
            for (i=0; i<numberToReceive; i++) {
                Atom &atom = system.addAtom();
                atom.position[0] = m_mpiReceiveBuffer[11*i + 0];
                atom.position[1] = m_mpiReceiveBuffer[11*i + 1];
                atom.position[2] = m_mpiReceiveBuffer[11*i + 2];

                atom.velocity[0] = m_mpiReceiveBuffer[11*i + 3];
                atom.velocity[1] = m_mpiReceiveBuffer[11*i + 4];
                atom.velocity[2] = m_mpiReceiveBuffer[11*i + 5];

                atom.initial_position[0] = m_mpiReceiveBuffer[11*i + 6];
                atom.initial_position[1] = m_mpiReceiveBuffer[11*i + 7];
                atom.initial_position[2] = m_mpiReceiveBuffer[11*i + 8];

                atom.setId(m_mpiReceiveBuffer[11*i + 9]);
                int atomicNumber = m_mpiReceiveBuffer[11*i + 10];
                atom.setType(AtomType::atomTypeFromAtomicNumber(atomicNumber));
            }
        }
    }
}
