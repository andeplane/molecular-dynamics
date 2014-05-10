#include "generator.h"
#include <atom.h>
#include <system.h>
#include <atomtype.h>

void Generator::generateFCC(System &system, double unitCellLength, vector<int> numberOfUnitCells, AtomType *atomType, double temperature)
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

