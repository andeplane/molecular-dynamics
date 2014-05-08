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
    int numberOfTimesteps = 500;
    double systemSize = 100;
    double delta = systemSize / 2;
    vector<double> systemLength(3,UnitConverter::lengthFromAngstroms(systemSize));

    simulator.initialize(0, vector<int>(3,1), systemLength);
    simulator.setTimestep(UnitConverter::timeFromSI(2e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
    cout << *potential << endl;
    FileManager fileManager;
//    fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());
    Atom &silicon = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Silicon));
    silicon.setPosition(UnitConverter::lengthFromAngstroms({1.825+delta, 1.826+delta, 1.826+delta}));

    Atom &oxygen1 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
    oxygen1.setPosition(UnitConverter::lengthFromAngstroms({0.896+delta,0.895+delta, 0.896+delta}));

    Atom &oxygen2 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
    oxygen2.setPosition(UnitConverter::lengthFromAngstroms({1.146+delta, 3.150+delta, 2.435+delta}));

    Atom &oxygen3 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
    oxygen3.setPosition(UnitConverter::lengthFromAngstroms({2.434+delta, 1.146+delta, 3.150+delta}));

    Atom &oxygen4 = simulator.system().addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
    oxygen4.setPosition(UnitConverter::lengthFromAngstroms({3.150+delta, 2.434+delta, 1.145+delta}));

    // simulator.system().atomManager().atoms().resetVelocityMaxwellian(UnitConverter::temperatureFromSI(0));

    simulator.system().atomManager().setGhostAtomsEnabled(false);
    cout << simulator.system() << endl;
    simulator.system().removeTotalMomentum();
    cout << simulator.system() << endl;
//    cout << "Before: " << endl << simulator.system().atomManager()<< endl;

//    simulator.step();

//    cout << "After: " << endl << simulator.system().atomManager()<< endl;
//    fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
//    fileManager.finalize();
//    return 0;

    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        simulator.step();


//        if(timestep%100 == 0) {
//            cout << "Timestep " << timestep;
//            cout << "   Total energy: " << simulator.sampler().calculateTotalEnergy(simulator.system());
//            cout << " (K=" << simulator.sampler().calculateKineticEnergy(simulator.system()) << "  P=" << simulator.sampler().calculatePotentialEnergy(simulator.system()) << ")" << endl;
//        } else {
//            cout << "Timestep " << timestep << endl;
//        }

    }
    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
    fileManager.finalize();
    return 0;
}

