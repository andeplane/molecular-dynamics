#include <vector>
#include <iostream>
#include <cmath>
using std::vector;
using std::cout; using std::endl;

#include <UnitTest++.h>
#include <cell.h>
#include <atom.h>
#include <atomtype.h>
#include <system.h>
#include <generator.h>
#include <simulator.h>
#include <integrators/velocityverlet.h>


SUITE(System) {
    TEST(AddRemoveAtoms) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        double cutOffDistance = 0.25;
        system.initialize(nodeIndex, numNodesVector, systemLength);

        CHECK_EQUAL(system.atomManager().atomCapacity(),1000);

        Atom *nullPointer = 0;
        Atom *atom = system.addAtom();
        CHECK_EQUAL(1, system.numberOfAtoms());
        system.addAtom();
        CHECK_EQUAL(2, system.numberOfAtoms());
        system.removeAtom(atom);
        CHECK_EQUAL(1, system.numberOfAtoms());
        system.removeAtom(atom);
        CHECK_EQUAL(0, system.numberOfAtoms());
        system.addAtom();
        CHECK_EQUAL(1, system.numberOfAtoms());
        system.addAtom();
        CHECK_EQUAL(2, system.numberOfAtoms());
        system.removeAtom(nullPointer);
        CHECK_EQUAL(2, system.numberOfAtoms());

        for(int i=0; i<1000; i++) {
            system.addAtom();
        }

        CHECK_EQUAL(2002, system.atomManager().atomCapacity());
        CHECK_EQUAL(1002, system.numberOfAtoms());

        system.removeAllAtoms();
        CHECK_EQUAL(0, system.numberOfAtoms());

        for(int i=0; i<1000; i++) {
            system.addAtom();
        }

        CHECK_EQUAL(1000, system.numberOfAtoms());
        CHECK_EQUAL(2002, system.atomManager().atomCapacity());
    }

    TEST(Cells) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        system.initialize(nodeIndex, numNodesVector, systemLength);

        Atom *atom = system.addAtom();
        CellData &cellData = system.atomManager().cellData();
        system.atomManager().setCutoffDistance(0.25);
        cellData = system.atomManager().cellData();
        Cell *cell = Cell::cellContainingAtom(atom,cellData);
        int cellIndex = Cell::cellIndexFromIJK(1,1,1,cellData);
        CHECK_EQUAL(cellIndex,cell->cellIndex());
    }

    TEST(IntegratorVelocityVerlet) {
        double timestep = 0.01;

        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        system.initialize(nodeIndex, numNodesVector, systemLength);

        vector<int> numberOfUnitCells(3,1);
        double fccLatticeConstant = 1.0;
        Atom *atom = system.addAtom();
        atom->setType(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        atom->velocity[0] = -1;
        atom->velocity[1] = 2;
        atom->velocity[2] = 3;
        vector<double> expectedPosition(3,0);

        VelocityVerlet integrator;
        for(int i=0; i<1000; i++) {
            expectedPosition[0] = fmod(atom->position[0] + atom->velocity[0]*timestep + systemLength[0],systemLength[0]);
            expectedPosition[1] = fmod(atom->position[1] + atom->velocity[1]*timestep + systemLength[1],systemLength[1]);
            expectedPosition[2] = fmod(atom->position[2] + atom->velocity[2]*timestep + systemLength[2],systemLength[2]);
            integrator.integrate(system,timestep);
            CHECK_ARRAY_CLOSE(expectedPosition,atom->position,3,1e-7);
        }

        atom->position[0] = 1; atom->position[1] = 1; atom->position[2] = 1;
        atom->velocity[0] = 1; atom->velocity[1] = 1; atom->velocity[2] = 1;
        expectedPosition[0] = fmod(atom->position[0] + atom->velocity[0]*timestep + systemLength[0],systemLength[0]);
        expectedPosition[1] = fmod(atom->position[1] + atom->velocity[1]*timestep + systemLength[1],systemLength[1]);
        expectedPosition[2] = fmod(atom->position[2] + atom->velocity[2]*timestep + systemLength[2],systemLength[2]);

        integrator.integrate(system,timestep);
        CHECK_ARRAY_CLOSE(expectedPosition,atom->position,3,1e-7);
    }
}
