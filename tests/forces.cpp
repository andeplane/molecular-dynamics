#include <vector>
using std::vector;


#include <UnitTest++.h>
#include <potentials/lennardjonespotential.h>
#include <includes.h>

SUITE(Forces) {
    TEST(ArgonForceTest) {
        Atom atom1(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        Atom atom2(AtomType::atomTypeFromAtomType(AtomTypes::Argon));

        atom1.setPosition(1, 0, 0);

        LennardJonesPotential potential;
        potential.twoParticleAction(&atom1,&atom2);
        CHECK_EQUAL(atom1.force[0], -atom2.force[0]);
        CHECK_EQUAL(atom1.force[1], -atom2.force[1]);
        CHECK_EQUAL(atom1.force[2], -atom2.force[2]);
    }
}
