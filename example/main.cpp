#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
#include <filemanager.h>
#include <statisticssampler.h>

using namespace std;

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    int numberOfTimesteps = 50;
    vector<double> systemLength(3,1.0);

    simulator.initialize(0, vector<int>(3,1), systemLength);
    simulator.setTimestep(0.02);

    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
    double latticeConstant = 1.54478708;
    Generator::generateFCC(simulator.system(),latticeConstant,vector<int>(3,5), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    simulator.system().removeTotalMomentum();

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

