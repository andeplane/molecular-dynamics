#include <iostream>
#include <simulator.h>
#include <generator.h>
#include <atomtype.h>
#include <vector>

using namespace std;

int main()
{
    double fccB = 1.54478708;
    vector<int> numNodesVector(3,1);
    vector<int> numUnitCells(3,1);
    vector<double> systemLength(3,1);
    numNodesVector[0] = 1; numNodesVector[1] = 1; numNodesVector[2] = 1;
    numUnitCells[0] = 1; numUnitCells[1] = 1; numUnitCells[2] = 1;

    systemLength[0] = numUnitCells[0]*fccB;
    systemLength[1] = numUnitCells[1]*fccB;
    systemLength[2] = numUnitCells[2]*fccB;

    Simulator simulator(0,numNodesVector, systemLength, 2.5);
    // simulator.system().addPotential(PotentialType::LennardJones);

    Generator generator;
    AtomType *type = AtomType::atomTypeFromAtomType(AtomTypes::Argon);

    // AtomType *type = AtomType::atomTypeFromAtomicNumber(18);
    // AtomType *type = AtomType::atomTypeFromAtomicNumber(9);

    generator.generateFCC(simulator.system(), fccB, numUnitCells, type);

    for(int i=0; i<10; i++) {
        simulator.step();
    }


    cout << "Finished" << endl;

    return 0;
}

