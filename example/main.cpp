#include <iostream>
#include <includes.h>
#include <random.h>
#include <atomiterators/atomiteratordefault.h>
#include <filemanager/filemanager.h>
#include <statisticssampler.h>
#include <unitconverter.h>
#include <statistics/temperaturesampler.h>
#include <statistics/kineticenergysampler.h>
#include <statistics/potentialenergysampler.h>
#include <statistics/totalenergysampler.h>
using namespace std;

void calculateDensity(System &system) {
    double totalMass = system.numberOfAtoms()*2.0/3*AtomType::atomTypeFromAtomicNumber(+AtomTypes::Oxygen)->mass() + system.numberOfAtoms()*1.0/3*AtomType::atomTypeFromAtomicNumber(+AtomTypes::Silicon)->mass();
    double massInGrams = UnitConverter::massToSI(totalMass) * 1000;
    vector<double> systemLength = system.systemLength();

    double totalVolume = systemLength[0]*systemLength[1]*systemLength[2];
    double volumeConversion = UnitConverter::lengthToSI(1.0)*UnitConverter::lengthToSI(1.0)*UnitConverter::lengthToSI(1.0);

    double totalVolumeSI = totalVolume*volumeConversion;
    double totalVolumeCC = totalVolumeSI*100*100*100;

    cout << "Mass in grams: " << massInGrams << endl;
    cout << "Total volume cc: " << totalVolumeCC << endl;
    cout << "Density [g/cc]: " << massInGrams / totalVolumeCC << endl;
}

int main()
{
    Random::setSeed(1);
    Simulator simulator;
    int numberOfTimesteps = 1000;

    simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({50, 50, 50}));
    simulator.setTimestep(UnitConverter::timeFromSI(2.0e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system()->addPotential(PotentialType::USCSilica);
//    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
//    Generator::generateFCC(simulator.system(), UnitConverter::lengthFromAngstroms(5.26), {10,10,10});

    FileManager fileManager;
    // fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());
    // Generator::addSiO4Molecule(simulator.system(), {25, 25, 25});
    Generator::generateBetaCrystabolite(simulator.system(),{5,5,5},UnitConverter::temperatureFromSI(300));
    shared_ptr<System> system = simulator.system();

    auto list = system->getChildrenOfClass<KineticEnergySampler>(true);
    cout << list.size() << endl;

    system->removeTotalMomentum();
    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        fileManager.saveMovieFrame(simulator.system()->atomManager().atoms().atoms(),simulator.system()->topology());
        simulator.step();
        if(timestep % 100 == 0) {
//            double energy = simulator.sampler().calculateTotalEnergy(simulator.system());
//            // double energy = simulator.sampler().calculatePotentialEnergy(simulator.system());
//            double energyEv = UnitConverter::energyToEv(energy);
//            double energyEvPerParticle = energyEv / simulator.system().numberOfAtoms();

//            double temperature = simulator.sampler().calculateTemperature(simulator.system());
//            double temperatureSI = UnitConverter::temperatureToSI(temperature);

            // cout << timestep << ": E=" << energyEvPerParticle << " eV, T=" << temperatureSI << endl;
            cout << timestep << endl;
        }
    }

    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    return 0;
}

