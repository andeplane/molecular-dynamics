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
    int numberOfTimesteps = 5;

    simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({100, 100, 100}));
    simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
    cout << *potential << endl;
    FileManager fileManager;
    fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());

    simulator.system().removeTotalMomentum();

    simulator.system().atomManager().setCutoffDistance(UnitConverter::lengthFromAngstroms(5.2));
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
    cout << "Momentum at end: " << simulator.sampler().calculateTotalMomentum(simulator.system()) << endl;
    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;

    fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
    fileManager.finalize();
    return 0;
}

