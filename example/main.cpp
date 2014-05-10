#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
#include <filemanager/filemanager.h>
#include <statisticssampler.h>
#include <unitconverter.h>
using namespace std;

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    int numberOfTimesteps = 50;

    simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({10, 10, 10}));
    simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
    cout << *potential << endl;
    FileManager fileManager;
    // fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());

    simulator.system().removeTotalMomentum();

    Generator::addSiO4Molecule(simulator.system(), {5,5,5});
    simulator.system().atomManager().setCutoffDistance(INFINITY);
    cout << simulator.system().atomManager() << endl;
    exit(1);
//    simulator.system().atomManager().atoms().iterate([](Atom &atom) {
//        atom.addVelocity(-5e-4,2e-3,1e-4);
//    });
    cout << "Before adding: " << simulator.system().atomManager().atoms() << endl << endl;
    vector<double> systemLength = simulator.system().systemLength();
    simulator.system().atomManager().atoms().iterate([&](Atom &atom) {
        for(int i=-1; i<=1; i++) {
            for(int j=-1; j<=1; j++) {
                for(int k=-1; k<=1; k++) {
                    if(i == 0 && j == 0 && k == 0) continue;
                    double x = atom.position[0] + i*systemLength[0];
                    double y = atom.position[1] + j*systemLength[1];
                    double z = atom.position[2] + k*systemLength[2];
                    simulator.system().addAtom(atom.type(), {x,y,z});
                }
            }
        }
    });

    cout << "After adding: " << simulator.system().atomManager().atoms() << endl << endl;

    exit(1);

    cout << simulator.system() << endl;
    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());

        if(timestep % 10 == 0) {
            // cout << timestep << "..";
            cout << timestep << " momentum: " << simulator.sampler().calculateTotalMomentum(simulator.system()) << endl;
        }

        fflush(stdout);

        simulator.step();
    }
    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    cout << simulator.system().atomManager().atoms() << endl;
    fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
    fileManager.finalize();
    return 0;
}

