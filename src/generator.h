#pragma once
#include <vector>
#include <unitconverter.h>
using std::vector;

class System; class AtomType;

class Generator
{
public:
    static void generateFCC(System &system, double latticeConstant, vector<int> numberOfUnitCells, AtomType *atomType, double temperature = 1.0);
    static void addSiO4Molecule(System &system, vector<double> SiPosition = {5,5,5}, double temperature = UnitConverter::temperatureFromSI(300), double tetrahedronScaling = 3.63);
};
