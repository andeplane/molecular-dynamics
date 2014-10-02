#pragma once
#include <vector>
#include <fstream>
using std::vector; using std::ifstream;

#include <utils/vec3.h>
#include <particles/atom.h>
class System;

class Mts0IO
{
public:
    Mts0IO(vector<int> numberOfCPUs);
    ~Mts0IO();
    void loadMts0(string mts0Directory, System &system);
    vector<int> numberOfCPUs() const;
    void setNumberOfCPUs(const vector<int> &numberOfCPUs);

private:
    void readData(ifstream *file, void *value);
    void readMts(char *filename, System &system);
    CompPhys::vec3 m_systemLength;
    vector<int> m_atomicNumberFromMts0AtomTypeMap;
    vector<double> m_nodeVectorIndex;
    vector<double> m_nodeOrigin;
    vector<double> m_nodeOffset;
    vector<int> m_numberOfCPUs;
};
