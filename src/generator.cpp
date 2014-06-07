#include "generator.h"
#include <atom.h>
#include <system.h>

void Generator::generateFCC(System &system, double unitCellLength, vector<int> numberOfUnitCells, shared_ptr<AtomType> atomType, double temperature)
{
    system.removeAllAtoms();

    vector<double> systemLength(3,0);
    systemLength[0] = numberOfUnitCells[0]*unitCellLength;
    systemLength[1] = numberOfUnitCells[1]*unitCellLength;
    systemLength[2] = numberOfUnitCells[2]*unitCellLength;
    system.setSystemLength(systemLength);

    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
    int atomID = 0;
    for(int x = 0; x < numberOfUnitCells[0]; x++) {
        for(int y = 0; y < numberOfUnitCells[1]; y++) {
            for(int z = 0; z < numberOfUnitCells[2]; z++) {
                for(int k = 0; k < 4; k++) {
                    Atom &atom = system.addAtom(atomType);
                    atom.position[0] = (x+xCell[k]) * unitCellLength;
                    atom.position[1] = (y+yCell[k]) * unitCellLength;
                    atom.position[2] = (z+zCell[k]) * unitCellLength;
                    atom.setId(atomID++);
                }
            }
        }
    }

    system.atomManager().atoms().resetVelocityMaxwellian(temperature);
}

void Generator::addSiO4Molecule(System &system, vector<double> SiPosition, double temperature, double tetrahedronScaling)
{
    double deltaX = SiPosition.at(0);
    double deltaY = SiPosition.at(1);
    double deltaZ = SiPosition.at(2);
    double D = tetrahedronScaling;
    system.addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Silicon),{deltaX, deltaY, deltaZ});
    system.addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen), {D/sqrt(3)  + deltaX, D/sqrt(3)  + deltaY, D/sqrt(3)  + deltaZ});
    system.addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen), {-D/sqrt(3) + deltaX, -D/sqrt(3) + deltaY, D/sqrt(3)  + deltaZ});
    system.addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen), {D/sqrt(3)  + deltaX, -D/sqrt(3) + deltaY, -D/sqrt(3) + deltaZ});
    system.addAtom(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen), {-D/sqrt(3) + deltaX, D/sqrt(3)  + deltaY, -D/sqrt(3) + deltaZ});
    system.atomManager().applyPeriodicBoundaryConditions(system.systemLength());

    system.atomManager().atoms().resetVelocityMaxwellian(temperature);
}

void Generator::generateBetaCrystabolite(System &system, vector<int> numberOfUnitCells, double temperature) {
    /*----------------------------------------------------------------------c
      Data for real beta-crystobalite SiO2
        NPU=24, NUA=8, and NUB=16
        Cubic unit cell A0 = B0 = C0 = 7.16 A = 13.5304 a.u.
        (R.W.G. Wyckoff, Crystal Structures, Vol. 1, P.318)
    c----------------------------------------------------------------------*/
    const int numberOfAtomsPerUnitCell = 24;
    const int numberOfSiliconPerUnitCell = 8;
    const int numberOfOxygenPerUnitCell = 16;

    double a0[3]={13.5304, 13.5304, 13.5304}; // In atomic units
    double pix[numberOfAtomsPerUnitCell][3]={
         {.255,.255,.255},{.755,.245,.745},{.245,.745,.755},{.745,.755,.245},
         {.992,.992,.992},{.492,.508,.008},{.508,.008,.492},{.008,.492,.508},
         {.125,.125,.125},{.625,.375,.875},{.375,.875,.625},{.875,.625,.375},
         {.660,.660,.060},{.160,.840,.940},{.160,.440,.340},{.560,.840,.340},
         {.660,.060,.660},{.340,.160,.440},{.340,.560,.840},{.940,.160,.840},
         {.060,.660,.660},{.840,.340,.560},{.840,.940,.160},{.440,.340,.160}};

    for (int j=0; j<numberOfAtomsPerUnitCell; j++) {
        bool atomIsSilicon = j < numberOfSiliconPerUnitCell;
        bool atomIsOxygen = !atomIsSilicon;

        shared_ptr<AtomType> atomType = AtomType::atomTypeFromAtomType(AtomTypes::Silicon);
        if(atomIsOxygen) atomType = AtomType::atomTypeFromAtomType(AtomTypes::Oxygen);

        for (int nZ=0; nZ<numberOfUnitCells[2]; nZ++) {
            for (int nY=0; nY<numberOfUnitCells[1]; nY++) {
                for (int nX=0; nX<numberOfUnitCells[0]; nX++) {
                    double x = (nX+pix[j][0])*a0[0];
                    double y = (nY+pix[j][1])*a0[1];
                    double z = (nZ+pix[j][2])*a0[2];
                    Atom &atom = system.addAtom(atomType, {x,y,z});
                }
            }
        }
    }

    vector<double> systemLength = {numberOfUnitCells[0]*a0[0], numberOfUnitCells[1]*a0[1], numberOfUnitCells[2]*a0[2]};
    system.setSystemLength(systemLength);

    system.atomManager().atoms().resetVelocityMaxwellian(temperature);
}
