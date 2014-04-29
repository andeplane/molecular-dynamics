#ifndef GENERATOR_H
#define GENERATOR_H
#include <vector>
using std::vector;

class System; class AtomType;

class Generator
{
public:
    Generator();
    void generateFCC(System &system, double latticeConstant, vector<int> numberOfUnitCells, AtomType *atomType);
};

#endif // GENERATOR_H
