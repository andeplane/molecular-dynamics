#pragma once
#include <vector>
#include <memory>
#include <utils/vec3.h>

using std::shared_ptr; using std::vector;

class System; class Atom; class AtomManager;

class Topology
{
private:
    int m_processorIndex;                                // MPI index
    int m_numProcessorsVector[3];                        // Number of processors in x,y,z direction
    int m_numProcessors;                                 // Total number of processors
    int m_neighborProcessorsIndices[6];                  // MPI index of the 6 face neighbors
    int m_parity[3];                                // Parity of each dimension (used to communicate in correct order over MPI)
    double m_shiftVector[6][3];                     // Contains the displacement in the physical space from this processor to each of the 6 face neighbors
    vector<int> m_processorCoordinates;                      // Three dimensional processor coordinate
    vector<double> m_processorLength;                    // Physical size of this processor
    CompPhys::vec3 m_systemLength;                  // System length in md units
    vector<double> m_origo;
    vector<vector<unsigned long> > m_moveQueue;     // Queue for atom indices to be moved to other processorr
    vector<double> m_mpiSendBuffer;
    vector<double> m_mpiReceiveBuffer;
    bool m_isInitialized;

    bool atomShouldBeCopied(Atom &atom, int &dimension, bool higher, double &cutoffDistance);
    bool atomDidChangeNode(Atom &atom, int &dimension, bool higher);
public:
    Topology();
    ~Topology();
    void initialize(int processorIndex, vector<int> numProcessorsVector, CompPhys::vec3 systemLength);
    void MPIMove(AtomManager &atomManager);         // Will move atoms between processors
    void copyGhostAtomsWithMPI(AtomManager &atomManager);
    int numProcessors() const;
    int processorIndex() const;
    CompPhys::vec3 systemLength() const;
    void setSystemLength(CompPhys::vec3 systemLength);
    vector<double> processorLength() const;
    vector<int> processorCoordinates() const;
    double maxSystemLength() const;
    vector<double> origo() const;
};
