#include "filemanager.h"
#include <topology.h>
#include <atom.h>
#include <atomtype.h>

FileManager::FileManager() :
    m_movieFile(NULL)
{

}

FileManager::~FileManager()
{
    if(isMovieFileOpen()) m_movieFile->close();
}

bool FileManager::isMovieFileOpen() const {
    return (m_movieFile!=NULL);
}

ofstream *FileManager::getMovieFile() const
{
    return m_movieFile;
}

void FileManager::finalize()
{
    if(isMovieFileOpen()) m_movieFile->close();
}

void FileManager::saveMovieFrame(const vector<Atom> &atoms, Topology &topology) {
    if(!isMovieFileOpen()) {
        char *filename = new char[10000];
        sprintf(filename,"movie_files/movie%04d.bin",topology.nodeIndex());
        m_movieFile = new ofstream(filename,std::ios::out | std::ios::binary);

        delete filename;
    }
    m_dataArray.clear();
    for(const Atom &atom : atoms) {
        double x = atom.position[0] + topology.origo()[0];
        double y = atom.position[1] + topology.origo()[1];
        double z = atom.position[2] + topology.origo()[2];
        double atomicNumber = atom.type()->atomicNumber();
        double atomID = atom.id();
        m_dataArray.push_back(x);
        m_dataArray.push_back(y);
        m_dataArray.push_back(z);
        m_dataArray.push_back(atomicNumber);
        m_dataArray.push_back(atomID);
    }

    unsigned long numberOfAtoms = atoms.size();

    m_movieFile->write (reinterpret_cast<char*>(&numberOfAtoms), sizeof(unsigned long));
    m_movieFile->write (reinterpret_cast<char*>(&m_dataArray.front()), 5*numberOfAtoms*sizeof(double));
}

void FileManager::loadState(string stateFolder, vector<Atom> &atoms) {

}

void FileManager::saveState(string stateFolder, vector<Atom> &atoms) {

}
