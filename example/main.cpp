#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
using namespace std;

int main()
{
    Random::setSeed(2);
    Simulator simulator;
    vector<double> systemLength(3,10);
    simulator.initialize(0, vector<int>(3,1), systemLength);
    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);

    // Generator::generateFCC(simulator.system(), 1, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    Atom &atom1 = simulator.system().addAtom();
    Atom &atom2 = simulator.system().addAtom();
    atom1.setType(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    atom2.setType(AtomType::atomTypeFromAtomType(AtomTypes::Argon));

    atom1.setPosition(2,0,0);

    cout << simulator.system().atomManager().atoms() << endl;
    simulator.step();
    cout << simulator.system().atomManager().atoms() << endl;

    return 0;
}

