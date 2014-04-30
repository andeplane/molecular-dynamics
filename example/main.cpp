#include <iostream>
#include <includes.h>

using namespace std;

int main()
{
    Simulator simulator;
    vector<double> systemLength(3,1);
    simulator.initialize(0, vector<int>(3,1), systemLength);
    Generator::generateFCC(simulator.system(), 1, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));

    simulator.system().atomManager().atoms().iterate([](Atom &atom, const int &atomIndex) {
        cout << "Atom with index " << atomIndex << ": " << atom << endl;
    });

    cout << "Num before moved: " << simulator.system().atomManager().numberOfAtoms() << endl;

    return 0;
}

