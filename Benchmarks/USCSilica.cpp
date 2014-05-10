#include <vector>
#include <iostream>
using std::vector;
using std::cout; using std::endl;

#include <UnitTest++.h>
#include <cell.h>
#include <atom.h>
#include <atomtype.h>
#include <system.h>
#include <generator.h>
#include <random.h>
#include <includes.h>
#include <filemanager/filemanager.h>

SUITE(USCSilica) {
    TEST(F77STATE) {
//        Random::setSeed(1);
//        Simulator simulator;
//        int numberOfTimesteps = 50;

//        simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({100, 100, 100}));
//        simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

//        USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
//        FileManager fileManager;
//        fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());

//        for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
//            if(timestep % 10 == 0) cout << timestep << "..";
//            fflush(stdout);

//            simulator.step();
//        }

        Random::setSeed(1);
        Simulator simulator;
        int numberOfTimesteps = 200;

        simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({100, 100, 100}));
        simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

        USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);

        FileManager fileManager;
        fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());

        simulator.system().removeTotalMomentum();

        for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
            // fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
            if(timestep % 10 == 0) cout << timestep << endl;

            simulator.step();
        }
        cout << "Momentum at end: " << simulator.sampler().calculateTotalMomentum(simulator.system()) << endl;
        cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    }
}
