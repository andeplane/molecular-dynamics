#include <vector>
#include <iostream>
#include <cmath>
using std::vector;
using std::cout; using std::endl;

#include <UnitTest++.h>
#include <includes.h>

SUITE(System) {
    TEST(AddRemoveAtoms) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        double cutOffDistance = 0.25;
        system.initialize(nodeIndex, numNodesVector, systemLength);

        CHECK_EQUAL(0, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());
        Generator::generateFCC(system,1.0, vector<int>(3,1),AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        CHECK_EQUAL(4, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());

        system.atomManager().atoms().iterate([&](Atom &atom, const int &atomIndex) {
            if(atomIndex % 2) {
                atom.setRemoved(true);
            }
        });

        CHECK_EQUAL(2, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());

        system.atomManager().addAtom();
        CHECK_EQUAL(3, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());

        system.atomManager().addGhostAtom();
        CHECK_EQUAL(3, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(1, system.atomManager().numberOfGhostAtoms());

        system.atomManager().removeGhostAtoms();
        CHECK_EQUAL(3, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());

        system.atomManager().addGhostAtom();
        CHECK_EQUAL(3, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(1, system.atomManager().numberOfGhostAtoms());

        system.atomManager().removeAllAtoms();
        CHECK_EQUAL(0, system.atomManager().numberOfAtoms());
        CHECK_EQUAL(0, system.atomManager().numberOfGhostAtoms());

        Generator::generateFCC(system, 1,vector<int>(3,10), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        CHECK_EQUAL(4000, system.atomManager().numberOfAtoms());

        Generator::generateFCC(system, 1,vector<int>(3,30), AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        CHECK_EQUAL(108000, system.atomManager().numberOfAtoms());
    }

    TEST(Cells) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        system.initialize(nodeIndex, numNodesVector, systemLength);

        Atom &atom = system.addAtom();
        CellData &cellData = system.atomManager().cellData();
        system.atomManager().setCutoffDistance(0.25);
        cellData = system.atomManager().cellData();

        Cell *cell = Cell::cellContainingAtom(atom,cellData);
        int cellIndex = Cell::cellIndexFromIJK(1,1,1,cellData);
        CHECK_EQUAL(cellIndex,cell->cellIndex());

        atom.setPosition(-0.1, -0.1, -0.1);
        cellData = system.atomManager().cellData(); // Update cell list
        cell = Cell::cellContainingAtom(atom,cellData);
        cellIndex = Cell::cellIndexFromIJK(0,0,0,cellData);
        CHECK_EQUAL(cellIndex,cell->cellIndex());

        atom.setPosition(0.25, 0.25, 0);
        cellData = system.atomManager().cellData(); // Update cell list
        cell = Cell::cellContainingAtom(atom,cellData);
        cellIndex = Cell::cellIndexFromIJK(2,2,1,cellData);
        CHECK_EQUAL(cellIndex,cell->cellIndex());

        atom.setPosition(0.25, 0.25, -0.2);
        cellData = system.atomManager().cellData(); // Update cell list
        cell = Cell::cellContainingAtom(atom,cellData);
        cellIndex = Cell::cellIndexFromIJK(2,2,0,cellData);
        CHECK_EQUAL(cellIndex,cell->cellIndex());
    }

    TEST(IntegratorVelocityVerlet) {
        double timestep = 0.01;

        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        system.initialize(nodeIndex, numNodesVector, systemLength);

        Atom &atom = system.addAtom();
        atom.setType(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        atom.setVelocity(-1,2,3);
        vector<double> expectedPosition(3,0);

        // Test with constant velocity
        VelocityVerlet integrator;
        for(int i=0; i<10000; i++) {
            expectedPosition[0] = fmod(atom.position[0] + atom.velocity[0]*timestep + systemLength[0],systemLength[0]);
            expectedPosition[1] = fmod(atom.position[1] + atom.velocity[1]*timestep + systemLength[1],systemLength[1]);
            expectedPosition[2] = fmod(atom.position[2] + atom.velocity[2]*timestep + systemLength[2],systemLength[2]);
            integrator.integrate(system,timestep);
            CHECK_ARRAY_CLOSE(expectedPosition,atom.position,3,1e-7);
            CHECK_EQUAL(false, atom.ghost());
            CHECK_EQUAL(1,system.numberOfAtoms());
        }

        atom.setPosition(1,1,1);
        atom.setVelocity(1,1,1);
        expectedPosition[0] = fmod(atom.position[0] + atom.velocity[0]*timestep + systemLength[0],systemLength[0]);
        expectedPosition[1] = fmod(atom.position[1] + atom.velocity[1]*timestep + systemLength[1],systemLength[1]);
        expectedPosition[2] = fmod(atom.position[2] + atom.velocity[2]*timestep + systemLength[2],systemLength[2]);
        integrator.integrate(system,timestep);
        CHECK_ARRAY_CLOSE(expectedPosition,atom.position,3,1e-7);

        // Now test with random velocity
        atom.resetMaxwellianVelocity(1.0);
        for(int i=0; i<10000; i++) {
            expectedPosition[0] = fmod(atom.position[0] + atom.velocity[0]*timestep + systemLength[0],systemLength[0]);
            expectedPosition[1] = fmod(atom.position[1] + atom.velocity[1]*timestep + systemLength[1],systemLength[1]);
            expectedPosition[2] = fmod(atom.position[2] + atom.velocity[2]*timestep + systemLength[2],systemLength[2]);
            integrator.integrate(system,timestep);
            CHECK_ARRAY_CLOSE(expectedPosition,atom.position,3,1e-7);
            CHECK_EQUAL(false, atom.ghost());
            CHECK_EQUAL(1,system.numberOfAtoms());
        }

        // Test with random velocity each timestep
        for(int i=0; i<10000; i++) {
            atom.resetMaxwellianVelocity(1.0);
            expectedPosition[0] = fmod(atom.position[0] + atom.velocity[0]*timestep + systemLength[0],systemLength[0]);
            expectedPosition[1] = fmod(atom.position[1] + atom.velocity[1]*timestep + systemLength[1],systemLength[1]);
            expectedPosition[2] = fmod(atom.position[2] + atom.velocity[2]*timestep + systemLength[2],systemLength[2]);
            integrator.integrate(system,timestep);
            CHECK_ARRAY_CLOSE(expectedPosition,atom.position,3,1e-7);
            CHECK_EQUAL(false, atom.ghost());
            CHECK_EQUAL(1,system.numberOfAtoms());
        }
    }
}
