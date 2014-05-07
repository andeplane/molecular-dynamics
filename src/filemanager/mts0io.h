#pragma once
#include <vector>
#include <fstream>
using std::vector; using std::ifstream;

#include <atom.h>
class System;

class Mts0IO
{
public:
    Mts0IO();
    vector<Atom> loadMts0(string mts0Directory, vector<int> numberOfCPUs, System &system);
private:
    void readData(ifstream *file, void *value);
    vector<int> m_atomicNumberFromMts0AtomTypeMap;
    void readMts(char *filename, System &system);
    vector<double> m_nodeVectorIndex;
    vector<double> m_nodeOrigin;
    vector<double> m_nodeOffset;
};
