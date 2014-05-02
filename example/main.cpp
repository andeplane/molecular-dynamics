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
    Generator::generateFCC(simulator.system(), latticeConstant, vector<int>(3,1), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//    Atom &atom1 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//    Atom &atom2 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//    atom2.setPosition(1.0,0,0);
//    atom2.setVelocity(-1,0,0);

    FileManager fileManager;
    int timesteps = 0;
    for(int i=0; i<500; i++) {
        if(i%100 == 0) {
            cout << i << endl;
        }
        simulator.step();
        fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        timesteps++;
    }

    fileManager.finalize();

    cout << "Saved " << timesteps << " movie frames." << endl;

    return 0;
}

