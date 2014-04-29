#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <vector>
using std::vector;

class System;
class Atom;

class Topology
{
private:
    int m_nodeIndex;                      // MPI index
    int m_numNodesVector[3];              // Number of processors in x,y,z direction
    int m_numNodes;                       // Total number of processors
    int m_neighborNodesIndices[6];        // MPI index of the 6 face neighbors
    int m_parity[3];                      // Parity of each dimension (used to communicate in correct order over MPI)
    double m_shiftVector[6][3];              // Contains the displacement in the physical space from this processor to each of the 6 face neighbors
    vector<int> m_nodeIndices;                 // Three dimensional processor coordinate
    vector<double> m_nodeLength;               // Physical size of this processor
    vector<double> m_systemLength;             // System length in md units
    vector<vector<int> > m_moveQueue;     // Queue for atom indices to be moved to other processor
    double *m_mpiSendBuffer;              // Temp array to save MPI data
    double *m_mpiReceiveBuffer;              // Temp array to save MPI data
    bool m_isInitialized;

    bool atomShouldBeCopied(Atom *atom, int &dimension, int &higher, double &cutoffDistance);
    bool atomDidChangeNode(Atom *atom, int &dimension, int &higher);
public:
    Topology();
    void initialize(int nodeIndex, vector<int> numNodesVector, vector<double> systemLength);
    void MPIMove(System &system);         // Will move atoms between processors
    void MPICopy(System &system, double cutoffDistance);         // Will copy ghost atoms between processors
    int numNodes() const;
    int nodeIndex() const;
    vector<double> systemLength() const;
    vector<double> nodeLength() const;
    vector<int> nodeIndices() const;
};

#endif // TOPOLOGY_H
