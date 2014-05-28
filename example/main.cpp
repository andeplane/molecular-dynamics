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
    int numberOfTimesteps = 1000;

    simulator.initialize(0, vector<int>(3,1), UnitConverter::lengthFromAngstroms({50, 50, 50}));
    simulator.setTimestep(UnitConverter::timeFromSI(1e-15));

    USCSIO2Potential *potential = (USCSIO2Potential*)simulator.system().addPotential(PotentialType::USCSilica);
//    LennardJonesPotential *potential = (LennardJonesPotential*)simulator.system().addPotential(PotentialType::LennardJones);
//    Generator::generateFCC(simulator.system(), UnitConverter::lengthFromAngstroms(5.26), {10,10,10});

    FileManager fileManager;
    // fileManager.loadMts0("/projects/andershaf_nanoporous_sio2_compressed_pore/test/heat/initial-crystal/mts0",{1,1,1},simulator.system());
    // Generator::addSiO4Molecule(simulator.system(), {25, 25, 25});
    simulator.system().addAtom(AtomType::atomTypeFromAtomicNumber(8), UnitConverter::lengthFromAngstroms({20,0,0}));
    simulator.system().addAtom(AtomType::atomTypeFromAtomicNumber(8), UnitConverter::lengthFromAngstroms({22,0,0}));
    simulator.system().addAtom(AtomType::atomTypeFromAtomicNumber(8), UnitConverter::lengthFromAngstroms({6,0,0}));
    simulator.system().addAtom(AtomType::atomTypeFromAtomicNumber(8), UnitConverter::lengthFromAngstroms({9,0,0}));
    simulator.system().addAtom(AtomType::atomTypeFromAtomicNumber(8), UnitConverter::lengthFromAngstroms({12,0,0}));
    // Generator::generateBetaCrystabolite(simulator.system(),{5,5,5});
    // simulator.system().setSystemLength({100, 100, 100});

    simulator.system().atomManager().atoms().iterate([](Atom &atom) {
//        double x = atom.position[0];
//        double y = atom.position[1];
//        double z = atom.position[2];
//        atom.setPosition(x+10,y+10,z+10);
        // atom.addVelocity( {1e-3, 2.2e-3, 0} );
    });

    // simulator.system().removeTotalMomentum();

    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        fileManager.saveMovieFrame(simulator.system().atomManager().atoms().atoms(),simulator.system().topology());
        cout << timestep << ": energy: " << simulator.sampler().calculateTotalEnergy(simulator.system()) << endl;

        simulator.step();
    }

    cout << "Momentum at end: " << simulator.sampler().calculateTotalMomentum(simulator.system()) << endl;
    cout << "Successfully computed " << numberOfTimesteps << " timesteps." << endl;
    return 0;
}

