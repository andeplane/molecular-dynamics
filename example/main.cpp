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
    FileManager fileManager;

    int timesteps = 0;
    vector<double> systemLength(3,3.0);

    simulator.initialize(0, vector<int>(3,1), systemLength);
    simulator.setTimestep(0.02);

    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
    double latticeConstant = 1.54478708;
    Generator::generateFCC(simulator.system(),latticeConstant,vector<int>(3,5), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
    simulator.system().removeTotalMomentum();

    for(int i=0; i<200; i++) {
        simulator.step();

        if(i%10 == 0) {
            cout << "Timestep " << i << "   Total energy: " << simulator.sampler().calculateTotalEnergy(simulator.system());
            cout << " (K=" << simulator.sampler().calculateKineticEnergy(simulator.system()) << "  P=" << simulator.sampler().calculatePotentialEnergy(simulator.system()) << ")" << endl;
        }

        // fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        timesteps++;
    }

    fileManager.finalize();

    cout << "Saved " << timesteps << " movie frames." << endl;

    return 0;
}

