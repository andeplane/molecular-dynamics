#include <statistics/neighborlist.h>
#include <atommanager.h>
#include <system.h>

double NeighborList::maxDistance() const
{
    return m_maxDistance;
}

void NeighborList::setMaxDistance(double maxDistance)
{
    m_maxDistance = maxDistance;
}

void NeighborList::action()
{
    m_list.clear();
    int numberOfCellsX = m_system->systemLength().x() / m_maxDistance;
    int numberOfCellsY = m_system->systemLength().y() / m_maxDistance;
    int numberOfCellsZ = m_system->systemLength().z() / m_maxDistance;
    numberOfCellsX = std::max(numberOfCellsX, 1); numberOfCellsY = std::max(numberOfCellsY, 1); numberOfCellsZ = std::max(numberOfCellsZ, 1);

    double cellLengthX = m_system->systemLength().x() / numberOfCellsX;
    double cellLengthY = m_system->systemLength().y() / numberOfCellsY;
    double cellLengthZ = m_system->systemLength().z() / numberOfCellsZ;

    numberOfCellsX += 2; numberOfCellsY += 2; numberOfCellsZ += 2; // Extra cells for the ghost atoms

    vector<vector<vector<vector<shared_ptr<Atom>>>>> cells(numberOfCellsX, vector<vector<vector<shared_ptr<Atom>>>>(numberOfCellsY, vector<vector<shared_ptr<Atom>>>(numberOfCellsZ)));
    auto atomSorter = [&](Atom &atom) {
        int cellIndexX = atom.position[0] / cellLengthX - 1; // Subtract 1 because of the extra ghost atom cells
        int cellIndexY = atom.position[1] / cellLengthY - 1;
        int cellIndexZ = atom.position[2] / cellLengthZ - 1;
        shared_ptr<Atom> atomPointer(&atom);

        cells.at(cellIndexX).at(cellIndexY).at(cellIndexZ).push_back(atomPointer);
    };

    m_system->atomManager().atoms().iterate(atomSorter);
    m_system->atomManager().ghostAtoms().iterate(atomSorter);
    for(int cellIndexX = 1; cellIndexX<numberOfCellsX-1; cellIndexX++) {
        for(int cellIndexY = 1; cellIndexY<numberOfCellsY-1; cellIndexY++) {
            for(int cellIndexZ = 1; cellIndexZ<numberOfCellsZ-1; cellIndexZ++) {
                auto &atomsInCell1 = cells.at(cellIndexX).at(cellIndexY).at(cellIndexZ);

                for(int dx = -1; dx<= 1; dx++) {
                    for(int dy = -1; dy<= 1; dy++) {
                        for(int dz = -1; dz<= 1; dz++) {
                            auto &atomsInCell2 = cells.at(cellIndexX+dx).at(cellIndexY+dy).at(cellIndexZ+dz);
                            for(auto atom1 : atomsInCell1) {
                                auto &neighborsOfThisAtom = m_list[atom1->uniqueId()];
                                for(auto atom2 : atomsInCell2) {
                                    neighborsOfThisAtom.push_back(atom2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

unordered_map<atomUniqueId, vector<shared_ptr<Atom>>> NeighborList::list() const
{
    return m_list;
}

NeighborList::NeighborList(shared_ptr<System> system, double maxDistance) :
    m_system(system),
    m_maxDistance(maxDistance)
{

}
