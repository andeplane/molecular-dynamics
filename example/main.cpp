#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
#include <filemanager.h>

using namespace std;

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    vector<double> systemLength(3,10);

    simulator.initialize(0, vector<int>(3,1), systemLength);
    simulator.setTimestep(0.02);
    // simulator.system().atomManager().setGhostAtomsEnabled(false);

    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
    double latticeConstant = 1.54478708;
    latticeConstant = 1.0;
    Generator::generateFCC(simulator.system(), latticeConstant, vector<int>(3,2), AtomType::atomTypeFromAtomType(AtomTypes::Argon));

//    Atom &atom1 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//    Atom &atom2 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//    atom2.setPosition(1.0,0,0);
//    atom2.setVelocity(-0.05,0,0);

    FileManager fileManager;
    int timesteps = 0;
//    cout << "Before timesteps" << endl;
//    cout << simulator.system().atomManager() << endl;
    cout << simulator.system() << endl;
    for(int i=0; i<10; i++) {
        if(i%100 == 0) {
            cout << i << endl;
        }
        simulator.step();
//        cout << "After step " << i << endl;
//        cout << simulator.system().atomManager() << endl;

        // fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        timesteps++;
    }

//    cout << "After step " << timesteps-1 << endl;
//    cout << simulator.system().atomManager() << endl;

    fileManager.finalize();

    cout << "Saved " << timesteps << " movie frames." << endl;

    return 0;
}

