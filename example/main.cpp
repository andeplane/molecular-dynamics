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
    int numberOfTimesteps = 2000;

    simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({100, 100, 100}));
    simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);

    FileManager fileManager;
    // fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());
    Generator::addSiO4Molecule(simulator.system(), {20, 20, 20});

    simulator.system().removeTotalMomentum();

    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        cout << timestep << ": energy: " << simulator.sampler().calculateTotalEnergy(simulator.system()) << endl;

        simulator.step();
    }

    cout << "Momentum at end: " << simulator.sampler().calculateTotalMomentum(simulator.system()) << endl;
    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    return 0;
}

