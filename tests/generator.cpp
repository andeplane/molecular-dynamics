#include <vector>
#include <iostream>
using std::vector;
using std::cout; using std::endl;

#include <UnitTest++.h>
#include <cell.h>
#include <particles/atom.h>
#include <atomtype.h>
#include <system.h>
#include <generator.h>

SUITE(Generator) {
    TEST(Generator) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        system.initialize(nodeIndex, numNodesVector, systemLength);

        vector<int> numberOfUnitCells(3,1);
        double fccLatticeConstant = 1.0;
        Generator::generateFCC(system, fccLatticeConstant, numberOfUnitCells, AtomType::atomTypeFromAtomType(AtomTypes::Argon));

        CHECK_EQUAL(4,system.numberOfAtoms());
        systemLength[0] = 1; systemLength[1] = 1; systemLength[2] = 1;
        CHECK_ARRAY_CLOSE(systemLength, system.topology().systemLength(), 3, 1e-7);

        system = System();
        system.initialize(nodeIndex, numNodesVector, systemLength);
        numberOfUnitCells[0] = 10; numberOfUnitCells[1] = 10; numberOfUnitCells[2] = 10;
        Generator::generateFCC(system, fccLatticeConstant, numberOfUnitCells, AtomType::atomTypeFromAtomType(AtomTypes::Argon));

        CHECK_EQUAL(4000,system.numberOfAtoms());
        systemLength[0] = 10; systemLength[1] = 10; systemLength[2] = 10;
        CHECK_ARRAY_CLOSE(systemLength, system.topology().systemLength(), 3, 1e-7);
    }
}
