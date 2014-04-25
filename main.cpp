#include <iostream>
#include <simulator.h>
#include <generator.h>

using namespace std;

int main()
{
    double fccB = 1.54478708;
    int numNodesVector[3];
    int numUnitCells[3];
    double systemLength[3];
    numNodesVector[0] = 1; numNodesVector[1] = 1; numNodesVector[2] = 1;
    numUnitCells[0] = 4; numUnitCells[1] = 4; numUnitCells[2] = 4;

    systemLength[0] = numUnitCells[0]*fccB;
    systemLength[1] = numUnitCells[1]*fccB;
    systemLength[2] = numUnitCells[2]*fccB;

    Simulator simulator(0,numNodesVector, systemLength, 2.5);
    simulator.system().addPotential(PotentialType::LennardJones);

    Generator generator;
    generator.generateFCC(simulator.system(), fccB, numUnitCells);

    for(int i=0; i<10; i++) {
        simulator.step();
    }


    cout << "Finished" << endl;

    return 0;
}

