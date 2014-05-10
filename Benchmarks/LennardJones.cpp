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

SUITE(LennardJones) {
    TEST(FCC) {
//        Random::setSeed(1);
//        Simulator simulator;
//        int numberOfTimesteps = 1000;
//        vector<double> systemLength(3,1.0);

//        simulator.initialize(0, vector<int>(3,1), systemLength);
//        simulator.setTimestep(0.02);

//        LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
//        double latticeConstant = 1.54478708;
//        Generator::generateFCC(simulator.system(),latticeConstant,vector<int>(3,5), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//        simulator.system().removeTotalMomentum();

//        for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
//            simulator.step();

//            if(timestep%100 == 0) {
//                cout << "Timestep " << timestep;
//                cout << "   Total energy: " << simulator.sampler().calculateTotalEnergy(simulator.system());
//                cout << " (K=" << simulator.sampler().calculateKineticEnergy(simulator.system()) << "  P=" << simulator.sampler().calculatePotentialEnergy(simulator.system()) << ")" << endl;
//            }
//        }
    }
}
