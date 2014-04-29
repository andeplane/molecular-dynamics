#include <vector>
using std::vector;


#include <UnitTest++.h>
#include <potentials/lennardjonespotential.h>
#include <atom.h>
#include <atomtype.h>
#include <system.h>

SUITE(Forces) {
    TEST(ArgonForceTest) {
//        Atom atom1(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//        Atom atom2(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
//        atom1.position[0] = 0;
//        atom2.position[0] = 0;

//        System system;
//        int nodeIndex = 0;
//        vector<int> numNodesVector(3,1);
//        vector<double> systemLength(3,1);
//        double cutOffDistance = 1.0;
//        system.initialize(nodeIndex, numNodesVector, systemLength, cutOffDistance);
//        system.atoms().push_back(atom1);
//        system.atoms().push_back(atom2);
//        system.updateCells();

//        LennardJonesPotential potential;
//        potential.calculateForces(system);
//        CHECK_EQUAL(atom1.force[0],-atom2.force[0]);
//        CHECK_EQUAL(atom1.force[1],-atom2.force[1]);
//        CHECK_EQUAL(atom1.force[2],-atom2.force[2]);
    }
}
