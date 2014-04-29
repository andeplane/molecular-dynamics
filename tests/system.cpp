#include <vector>
#include <iostream>
using std::vector;
using std::cout; using std::endl;
#include <UnitTest++.h>
#include <atom.h>
#include <atomtype.h>
#include <system.h>

SUITE(System) {
    TEST(AddRemoveAtoms) {
        System system;
        int nodeIndex = 0;
        vector<int> numNodesVector(3,1);
        vector<double> systemLength(3,1);
        double cutOffDistance = 0.25;
        system.initialize(nodeIndex, numNodesVector, systemLength, cutOffDistance);

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
}
