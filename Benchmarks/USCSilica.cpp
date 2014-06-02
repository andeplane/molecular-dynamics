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
        int numberOfTimesteps = 1000;

        simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({100, 100, 100}));
        simulator.setTimestep(UnitConverter::timeFromSI(2e-15));

        USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
        Generator::generateBetaCrystabolite(simulator.system(), {5,5,5}, UnitConverter::temperatureFromSI(300));
        simulator.system().removeTotalMomentum();

        for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
            simulator.step();
            if(timestep % 100 == 0) {
                double energy = simulator.sampler().calculateTotalEnergy(simulator.system());
                // double energy = simulator.sampler().calculatePotentialEnergy(simulator.system());
                double energyEv = UnitConverter::energyToEv(energy);
                double energyEvPerParticle = energyEv / simulator.system().numberOfAtoms();

                double temperature = simulator.sampler().calculateTemperature(simulator.system());
                double temperatureSI = UnitConverter::temperatureToSI(temperature);

                cout << timestep << ": E=" << energyEvPerParticle << " eV, T=" << temperatureSI << endl;
            }
        }
        cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    }
}
