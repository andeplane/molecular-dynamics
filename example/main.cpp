#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
using namespace std;

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    vector<double> systemLength(3,1);

    simulator.initialize(0, vector<int>(3,1), systemLength);
    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
    double latticeConstant = 1.54478708;
    Generator::generateFCC(simulator.system(), latticeConstant, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    systemLength = simulator.system().systemLength();

    for(int i=0; i<1000; i++) {
//        if(i%100 == 0) {
//            cout << i << endl;
//        }
        simulator.step();
    }

    return 0;
}

