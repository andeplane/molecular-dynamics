#include "generator.h"
#include <atom.h>
#include <system.h>
#include <atomtype.h>

Generator::Generator()
{

}

void Generator::generateFCC(System &system, double unitCellLength, vector<int> numberOfUnitCells, AtomType *atomType)
{
    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
    int atomID = 0;
    for(int x = 0; x < numberOfUnitCells[0]; x++) {
        for(int y = 0; y < numberOfUnitCells[1]; y++) {
            for(int z = 0; z < numberOfUnitCells[2]; z++) {
                for(int k = 0; k < 4; k++) {
                    // Set positions and type
                    Atom *atom = system.addAtom();
                    atom->setType(atomType);
                    atom->position[0] = (x+xCell[k]) * unitCellLength;
                    atom->position[1] = (y+yCell[k]) * unitCellLength;
                    atom->position[2] = (z+zCell[k]) * unitCellLength;
                    atom->setId(atomID++);
                }
            }
        }
    }
}
