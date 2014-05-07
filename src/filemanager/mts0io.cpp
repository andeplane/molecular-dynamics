#include <filemanager/mts0io.h>
#include <system.h>
#include <atommanager.h>

Mts0IO::Mts0IO(vector<int> numberOfCPUs) :
    m_numberOfCPUs(numberOfCPUs)
{
    m_atomicNumberFromMts0AtomTypeMap.resize(8);
    m_atomicNumberFromMts0AtomTypeMap[1] = +AtomTypes::Silicon;
    m_atomicNumberFromMts0AtomTypeMap[2] = +AtomTypes::Oxygen;
    m_atomicNumberFromMts0AtomTypeMap[3] = +AtomTypes::Hydrogen;
    m_atomicNumberFromMts0AtomTypeMap[4] = +AtomTypes::Oxygen;
    m_atomicNumberFromMts0AtomTypeMap[5] = +AtomTypes::Sodium;
    m_atomicNumberFromMts0AtomTypeMap[6] = +AtomTypes::Chlorine;
    m_systemLength.resize(3,0);

    if(numberOfCPUs.size()<3) {
        cout << "Error, numberOfCPUs only contains " << numberOfCPUs.size() << " elements. It needs to know how many cpu's in all three dimensions." << endl;
    }
}

Mts0IO::~Mts0IO()
{
    m_systemLength.clear();
    m_nodeOffset.clear();
    m_nodeOrigin.clear();
    m_nodeVectorIndex.clear();
    m_atomicNumberFromMts0AtomTypeMap.clear();
    m_numberOfCPUs.clear();
}

void Mts0IO::readData(ifstream *file, void *value) {
    int N;
    file->read (reinterpret_cast<char*>(&N), sizeof(int));
    file->read (reinterpret_cast<char*>(value), N);
    file->read (reinterpret_cast<char*>(&N), sizeof(int));
}

void Mts0IO::readMts(char *filename, System &system) {
    ifstream *file = new ifstream();
    file->open(filename, std::ios::in | std::ios::binary);
    if (!*file) {
        std::cerr << "Error in Mts0IO::read_mts(): Failed to open file " << filename << std::endl;
        exit(1);
    }

    int numberOfAtoms;
    readData(file, &numberOfAtoms);

    double *phaseSpace = new double[6*numberOfAtoms];
    double *tmpAtomData = new double[numberOfAtoms];

    readData(file, tmpAtomData);
    readData(file, phaseSpace);
    vector<vector<vector<double> > > hMatrix;

    hMatrix.resize(2,vector<vector<double> >(3,{0,0,0}));

    double *tmpHMatrix = new double[18];
    readData(file,tmpHMatrix);
    int count = 0;
    for(int k=0;k<2;k++) {
        for(int j=0;j<3;j++) {
            for(int i=0;i<3;i++) {
                hMatrix.at(k).at(i).at(j) = float(tmpHMatrix[count++]);
            }
        }
    }

    m_systemLength.at(0) = hMatrix[0][0][0];;
    m_systemLength.at(1) = hMatrix[0][1][1];
    m_systemLength.at(2) = hMatrix[0][2][2];

    for(int i=0;i<numberOfAtoms;i++) {
        int atomType = int(tmpAtomData[i]);
        unsigned long atomId = (tmpAtomData[i]-atomType)*1e11 + 1e-5; // Handle roundoff errors from 2 -> 1.99999999 -> 1
        int atomicNumber = m_atomicNumberFromMts0AtomTypeMap.at(atomType);
        double x = (phaseSpace[3*i+0] + m_nodeOrigin[0])*m_systemLength[0]; // Add node origin since all positions are local on each cpu
        double y = (phaseSpace[3*i+1] + m_nodeOrigin[1])*m_systemLength[1]; // Also multiply by system length because all positions are between 0 and 1 [1 beeing the maximum coordinate if one processor]
        double z = (phaseSpace[3*i+2] + m_nodeOrigin[2])*m_systemLength[2];
        double vx = phaseSpace[3*(numberOfAtoms + i)+0];
        double vy = phaseSpace[3*(numberOfAtoms + i)+1];
        double vz = phaseSpace[3*(numberOfAtoms + i)+2];

        Atom &atom = system.addAtom(AtomType::atomTypeFromAtomicNumber(atomicNumber));
        atom.setPosition(x,y,z);
        atom.setVelocity(vx,vy,vz);
        atom.setId(atomId);
    }

    file->close();
    hMatrix.clear();
    delete tmpHMatrix;
    delete phaseSpace;
    delete tmpAtomData;
    delete file;
}

void Mts0IO::loadMts0(string mts0Directory, System &system)
{
    int nx = m_numberOfCPUs.at(0);
    int ny = m_numberOfCPUs.at(1);
    int nz = m_numberOfCPUs.at(2);

    int numCPUs = nx*ny*nz;
    m_nodeVectorIndex.resize(3,0);
    m_nodeOrigin.resize(3,0);
    m_nodeOffset.resize(3,0);

    char filename[10000];
    m_nodeOffset[0] = 1.0/nx;
    m_nodeOffset[1] = 1.0/ny;
    m_nodeOffset[2] = 1.0/nz;

    for(int cpuID=0; cpuID<numCPUs; cpuID++) {
        m_nodeVectorIndex[0] = cpuID/(ny*nz);   // Node id in x-direction
        m_nodeVectorIndex[1] = (cpuID/nz) % ny; // Node id in y-direction
        m_nodeVectorIndex[2] = cpuID % nz; 	  // Node id in z-direction

        m_nodeOrigin[0] = m_nodeOffset[0]*m_nodeVectorIndex[0]; // Displacement in x-direction
        m_nodeOrigin[1] = m_nodeOffset[1]*m_nodeVectorIndex[1]; // Displacement in y-direction
        m_nodeOrigin[2] = m_nodeOffset[2]*m_nodeVectorIndex[2]; // Displacement in z-direction

        sprintf(filename,"%s/mt%04d",mts0Directory.c_str(), cpuID);
        readMts(filename,system);
    }

    system.setSystemLength(m_systemLength);
}
