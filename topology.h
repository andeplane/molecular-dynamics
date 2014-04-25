#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <vector>
using std::vector;

class System;
class Atom;

class Topology
{
private:
    int m_nodeIndices[3];                 // Three dimensional processor coordinate
    int m_nodeIndex;                      // MPI index
    int m_numNodesVector[3];              // Number of processors in x,y,z direction
    int m_numNodes;                       // Total number of processors
    int m_neighborNodesIndices[6];        // MPI index of the 6 face neighbors
    int m_parity[3];                      // Parity of each dimension (used to communicate in correct order over MPI)
    double m_nodeLength[3];               // Physical size of this processor
    double m_shiftVector[6][3];              // Contains the displacement in the physical space from this processor to each of the 6 face neighbors
    double m_systemLength[3];             // System length in md units
    vector<vector<int> > m_moveQueue;     // Queue for atom indices to be moved to other processor
    double *m_mpiSendBuffer;              // Temp array to save MPI data
    double *m_mpiReceiveBuffer;              // Temp array to save MPI data

    bool atomShouldBeCopied(Atom &atom, int &dimension, int &higher, double &cutoffDistance);
    bool atomDidChangeNode(Atom &atom, int &dimension, int &higher);
public:
    void initialize(int nodeIndex, int numNodesVector[3], double systemLength[3]);
    void MPIMove(System &system);         // Will move atoms between processors
    void MPICopy(System &system, double cutoffDistance);         // Will copy ghost atoms between processors
};

#endif // TOPOLOGY_H
