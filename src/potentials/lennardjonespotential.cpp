#include <potentials/lennardjonespotential.h>
#include <system.h>

LennardJonesPotential::LennardJonesPotential() :
    m_sigma(1),
    m_epsilon(1)
{

}

double LennardJonesPotential::sigma() const
{
    return m_sigma;
}

void LennardJonesPotential::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJonesPotential::epsilon() const
{
    return m_epsilon;
}

void LennardJonesPotential::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}
#include <iostream>
using namespace std;
void LennardJonesPotential::calculateForces(AtomManager &atomList)
{
//    atomList.updateCells();
//    double dr[3];
//    double cutoffDistanceSquared = atomList.cutoffDistance()*atomList.cutoffDistance();
//    cout << "will do stuff" << endl;
//    for(int cellX=1; cellX<=SystemCell::numberOfCellsWithoutGhostCells[0]; cellX++) {
//        for(int cellY=1; cellY<=SystemCell::numberOfCellsWithoutGhostCells[1]; cellY++) {
//            for(int cellZ=1; cellZ<=SystemCell::numberOfCellsWithoutGhostCells[2]; cellZ++) {
//                SystemCell &cell1 = atomList.cells().at(SystemCell::cellIndexFromIJK(cellX, cellY, cellZ));

//                for(int cell2X=cellX-1; cell2X<=cellX+1; cell2X++) {
//                    for(int cell2Y=cellY-1; cell2Y<=cellY+1; cell2Y++) {
//                        for(int cell2Z=cellZ-1; cell2Z<=cellZ+1; cell2Z++) {
//                            SystemCell &cell2 = atomList.cells().at(SystemCell::cellIndexFromIJK(cell2X, cell2Y, cell2Z));

//                            for(Atom *atom1 : cell1.atoms()) {
//                                for(Atom *atom2: cell2.atoms()) {
//                                    if(atom1->id() < atom2->id()) continue; // Newton's 3rd law

//                                    dr[0] = atom1->position[0] - atom2->position[0];
//                                    dr[1] = atom1->position[1] - atom2->position[1];
//                                    dr[2] = atom1->position[2] - atom2->position[2];
//                                    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
//                                    if (dr2<cutoffDistanceSquared) {
//                                        double dr2Inverse = 1.0/dr2;
//                                        double dr6Inverse = dr2Inverse*dr2Inverse*dr2Inverse;
//                                        double force = (2*dr6Inverse-1)*dr6Inverse*dr2Inverse*24;
//                                        cout << force << endl;

//                                        atom1->force[0] += force*dr[0];
//                                        atom1->force[1] += force*dr[1];
//                                        atom1->force[2] += force*dr[2];
//                                        if(!atom2->isGhost()) {
//                                            atom2->force[0] -= force*dr[0];
//                                            atom2->force[1] -= force*dr[1];
//                                            atom2->force[2] -= force*dr[2];
//                                        }
//                                    }

//                                }
//                            } // Loop atoms

//                        }
//                    }
//                } // Loop neighbor cells
//            }
//        }
//    }
}
