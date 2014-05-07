#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
#include <filemanager/filemanager.h>
#include <statisticssampler.h>

using namespace std;

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    int numberOfTimesteps = 50;
    vector<double> systemLength(3,1.0);

    simulator.initialize(0, vector<int>(3,1), systemLength);
    simulator.setTimestep(1.0);

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
    FileManager fileManager;
    fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());
    fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
    fileManager.finalize();
    return 0;

    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        simulator.step();

//        if(timestep%100 == 0) {
//            cout << "Timestep " << timestep;
//            cout << "   Total energy: " << simulator.sampler().calculateTotalEnergy(simulator.system());
//            cout << " (K=" << simulator.sampler().calculateKineticEnergy(simulator.system()) << "  P=" << simulator.sampler().calculatePotentialEnergy(simulator.system()) << ")" << endl;
//        } else {
//            cout << "Timestep " << timestep << endl;
//        }

    }
    return 0;
}

