#pragma once
#include <vector>
using std::vector;

class System; class AtomType;

class Generator
{
public:
    static void generateFCC(System &system, double latticeConstant, vector<int> numberOfUnitCells, AtomType *atomType, double temperature = 1.0);
};
