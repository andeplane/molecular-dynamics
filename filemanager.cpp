#include "filemanager.h"
#include <topology.h>
#include <atom.h>

FileManager::FileManager() :
    m_movieFile(NULL)
{

}

bool FileManager::isMovieFileOpen() const {
    return (m_movieFile!=NULL);
}

ofstream *FileManager::getMovieFile() const
{
    return m_movieFile;
}

void FileManager::saveMovieFrame(vector<Atom> &atoms, Topology &topology) {
    if(!isMovieFileOpen()) {
        char *filename = new char[10000];
        sprintf(filename,"movie_files/movie%04d.bin",topology.nodeIndex());
        m_movieFile = new ofstream(filename,std::ios::out | std::ios::binary);

        delete filename;
    }

    m_dataArray.clear();

    for(Atom &atom : atoms) {
        double x = atom.position[0] + topology.systemLength()[0];
        double y = atom.position[1] + topology.systemLength()[1];
        double z = atom.position[2] + topology.systemLength()[2];
        double atomType = atom.type();
        double atomID = atom.id();
        m_dataArray.push_back(x);
        m_dataArray.push_back(y);
        m_dataArray.push_back(z);
        m_dataArray.push_back(atomType);
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
