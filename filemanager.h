#ifndef FILEMANAGER_H
#define FILEMANAGER_H
#include <string>
#include <fstream>
#include <vector>

using std::string;
using std::vector;
using std::ofstream;

class System; class Atom; class Topology;

class FileManager
{
private:
    ofstream *m_movieFile;
    vector<double> m_dataArray;
protected:
    bool isMovieFileOpen() const;
public:
    FileManager();
    void loadState(string stateFolder, vector<Atom> &atoms);
    void saveState(string stateFolder, vector<Atom> &atoms);
    ofstream *getMovieFile() const;
    void setMovieFile(FILE *value);

    void saveMovieFrame(vector<Atom> &atoms, Topology &topology);
};

#endif // FILEMANAGER_H
