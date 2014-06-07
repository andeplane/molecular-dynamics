#pragma once
#include <vector>
#include <unitconverter.h>
#include <atomtype.h>
using std::vector;

class System; class AtomType;

class Generator
{
public:
    static void generateFCC(System &system, double latticeConstant = UnitConverter::lengthFromAngstroms(5.26), vector<int> numberOfUnitCells = {10, 10, 10}, shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::Argon), double temperature = UnitConverter::temperatureFromSI(100));
    static void addSiO4Molecule(System &system, vector<double> SiPosition = {5,5,5}, double temperature = UnitConverter::temperatureFromSI(300), double tetrahedronScaling = 3.63);
    static void generateBetaCrystabolite(shared_ptr<System> system, vector<int> numberOfUnitCells, double temperature = UnitConverter::temperatureFromSI(300));
};
