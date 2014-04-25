#include "topology.h"

#include <system.h>
#include <atom.h>
#include <algorithm>
using std::swap;

void Topology::initialize(int nodeIndex, int numNodesVector[3], double systemLength[3])
{
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/
    m_mpiReceiveBuffer = new double[1000000];
    m_mpiSendBuffer = new double[1000000];
    m_nodeIndex = nodeIndex;

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

bool Topology::atomShouldBeCopied(Atom &atom, int &dimension, int &higher, double &cutoffDistance) {
  if (higher == 0) return atom.position[dimension] < cutoffDistance;
  else return atom.position[dimension] > m_nodeLength[dimension]-cutoffDistance;
}

bool Topology::atomDidChangeNode(Atom &atom, int &dimension, int &higher) {
    if (higher == 0) return atom.position[dimension] < 0.0;
    else return atom.position[dimension] > m_nodeLength[dimension];
}

void Topology::MPICopy(System &system, double cutoffDistance)
{
    int numNewGhostAtoms = 0;
    system.setFirstGhostAtomIndex(system.atoms().size()); // First ghost atom will be the next

    for(vector<int> &queue : m_moveQueue) {
        queue.clear();
    }

    for(int dimension=0;dimension<3;dimension++) {
        for(int atomIndex=0; atomIndex<system.atoms().size(); atomIndex++) {
            Atom &atom = system.atoms().at(atomIndex);

            for(int higher=0; higher<=1; higher++) {
                int localNodeID = 2*dimension + higher;
                if(atomShouldBeCopied(atom,dimension,higher,cutoffDistance)) m_moveQueue.at(localNodeID).push_back(atomIndex);
            }
        }

        /* Loop through higher and lower node in this dimension */
        for(int higher=0;higher<2;higher++) {
            int localNodeID = 2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();
            int numberToRecieve = numberToSend;

            int i = 0;
            for (int atomIndex : m_moveQueue.at(localNodeID)) {
                Atom &atom = system.atoms().at(atomIndex);

                /* Shift the coordinate origin */
                m_mpiSendBuffer[3*i+0] = atom.position[0]-m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[3*i+1] = atom.position[1]-m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[3*i+2] = atom.position[2]-m_shiftVector[localNodeID][2];
                i++;
            }

            memcpy(m_mpiReceiveBuffer,m_mpiSendBuffer,3*numberToRecieve*sizeof(double));

            for(int i=0; i<numberToRecieve; i++) {
                Atom ghostAtom;
                ghostAtom.position[0] = m_mpiReceiveBuffer[3*i+0];
                ghostAtom.position[1] = m_mpiReceiveBuffer[3*i+1];
                ghostAtom.position[2] = m_mpiReceiveBuffer[3*i+2];

                system.atoms().push_back(ghostAtom);
            }

            numNewGhostAtoms += numberToRecieve;
        }
    }
}

void Topology::MPIMove(System &system) {
    int numNewAtoms = 0;
    system.atoms().erase(system.atoms().begin()+system.firstGhostAtomIndex(),system.atoms().end()); // Delete all ghost atoms

    for(vector<int> &queue : m_moveQueue) {
        queue.clear();
    }

    for(int dimension=0;dimension<3;dimension++) {
        for(int atomIndex=0; atomIndex<system.atoms().size(); atomIndex++) {
            Atom &atom = system.atoms().at(atomIndex);
            int nodeLower = 2*dimension;
            int nodeHigher = 2*dimension+1;
            /* Register a to-be-copied atom in move_queue[kul|kuh][] */
            if (!atom.moved()) { /* Don't scan moved-out atoms */
                if(atomDidChangeNode(atom,dimension,nodeLower)) {
                    m_moveQueue.at(nodeLower).push_back(atomIndex);
                } else if (atomDidChangeNode(atom, dimension, nodeHigher)) {
                    m_moveQueue.at(nodeHigher).push_back(atomIndex);
                }
            }
        }

        for (int higher=0; higher<=1; higher++) {
            int localNodeID=2*dimension+higher;
            int numberToSend = m_moveQueue.at(localNodeID).size();

            int numberToReceive = numberToSend;

            /* Message buffering */
            int i = 0;
            for (int atomIndex : m_moveQueue.at(localNodeID)) {
                Atom &atom = system.atoms().at(atomIndex);

                /* Shift the coordinate origin */
                m_mpiSendBuffer[11*(i-1)    + 0] = atom.position[0] - m_shiftVector[localNodeID][0];
                m_mpiSendBuffer[11*(i-1)    + 1] = atom.position[1] - m_shiftVector[localNodeID][1];
                m_mpiSendBuffer[11*(i-1)    + 2] = atom.position[2] - m_shiftVector[localNodeID][2];
                m_mpiSendBuffer[11*(i-1)+ 3 + 0] = atom.velocity[0];
                m_mpiSendBuffer[11*(i-1)+ 3 + 1] = atom.velocity[1];
                m_mpiSendBuffer[11*(i-1)+ 3 + 2] = atom.velocity[2];
                m_mpiSendBuffer[11*(i-1)+ 6 + 0] = atom.initial_position[0];
                m_mpiSendBuffer[11*(i-1)+ 6 + 1] = atom.initial_position[1];
                m_mpiSendBuffer[11*(i-1)+ 6 + 2] = atom.initial_position[2];
                m_mpiSendBuffer[11*(i-1) + 9]    = (double)atom.id();
                m_mpiSendBuffer[11*(i-1) + 10]   = (double)atom.type();
                atom.setMoved(true);
                i++;
            }

            memcpy(m_mpiReceiveBuffer,m_mpiSendBuffer,11*numberToReceive*sizeof(double));

            /* Message storing */
            for (i=0; i<numberToReceive; i++) {
                Atom atom;
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
                atom.setType(m_mpiReceiveBuffer[11*i + 10]);
            }

            /* Increment the # of new atoms */
            numNewAtoms += numberToReceive;
        }
    }

    for(int atomIndex=0; atomIndex<system.atoms().size(); atomIndex++) {
        Atom &atom = system.atoms().at(atomIndex);
        if(atom.moved()) {
            Atom &lastAtom = system.atoms().back();
            swap(atom,lastAtom);
            system.atoms().pop_back();
            atomIndex--;
        }
    }
}
