#include <iostream>
#include <includes.h>

using namespace std;

int main()
{
    Simulator simulator;
    vector<double> systemLength(3,1);
    simulator.initialize(0, vector<int>(3,1), systemLength);
    Generator::generateFCC(simulator.system(), 1, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    cout << "Number of atoms: " << simulator.system().numberOfAtoms() << endl;
    return 0;
}

