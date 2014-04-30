#include <iostream>
#include <includes.h>
#include <random.h>

using namespace std;

int main()
{
    Random::setSeed(2);
    Simulator simulator;
    vector<double> systemLength(3,1);
    simulator.initialize(0, vector<int>(3,1), systemLength);
    Generator::generateFCC(simulator.system(), 1, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    cout << simulator.system().atomManager().atoms() << endl;
    simulator.step();
    cout << simulator.system().atomManager().atoms() << endl;
    return 0;
}

