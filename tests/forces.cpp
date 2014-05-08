#include <vector>
using std::vector;
#include <UnitTest++.h>
#include <potentials/potentials.h>
#include <includes.h>

void lennardJones(Atom &atom1, Atom &atom2, vector<double> &forceVector, double cutoffDistance) {
    forceVector.resize(3,0);
    forceVector[0] = 0; forceVector[1] = 0; forceVector[2] = 0;
    double dr[3];
    dr[0] = atom1.position[0]-atom2.position[0];
    dr[1] = atom1.position[1]-atom2.position[1];
    dr[2] = atom1.position[2]-atom2.position[2];

    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    double cutoffDistanceSquared = cutoffDistance*cutoffDistance;
    if (dr2<cutoffDistanceSquared) {
        double dr2_inverse = 1.0/dr2;
        double dr6_inverse = dr2_inverse*dr2_inverse*dr2_inverse;
        double force = (2*dr6_inverse-1)*dr6_inverse*dr2_inverse*24;
        forceVector[0] = force*dr[0];
        forceVector[1] = force*dr[1];
        forceVector[2] = force*dr[2];
    }
}

void setPositionsAndResetForces(Atom &atom1, Atom &atom2, double r1x, double r1y, double r1z, double r2x, double r2y, double r2z) {
    atom1.resetForce(); atom2.resetForce();
    atom1.setPosition(r1x,r1y,r1z);
    atom2.setPosition(r2x,r2y,r2z);
}

SUITE(Forces) {
    TEST(ArgonForceTest) {
        Atom atom1(AtomType::atomTypeFromAtomType(AtomTypes::Argon));
        Atom atom2(AtomType::atomTypeFromAtomType(AtomTypes::Argon));

        LennardJonesPotential potential;
        potential.setSigma(1.0);
        potential.setEpsilon(1.0);

        atom1.setPosition(0.78432233, -0.112332, 1.22315322);
        potential.twoParticleAction(&atom1,&atom2);
        CHECK_EQUAL(atom1.force[0], -atom2.force[0]);
        CHECK_EQUAL(atom1.force[1], -atom2.force[1]);
        CHECK_EQUAL(atom1.force[2], -atom2.force[2]);
        vector<double> forceVector(3,0);

        setPositionsAndResetForces(atom1,atom2, 1,0,0, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, -1,0,0, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 0,1,0, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 0,-1,0, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 0,0,1, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 0,0,1, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 10,10,10, 0,0,0);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        setPositionsAndResetForces(atom1,atom2, 10,10,10, 5,4,3);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);

        potential.setCutoffDistance(INFINITY);
        setPositionsAndResetForces(atom1,atom2, -10,-10,10, -5,4,3);
        lennardJones(atom1,atom2,forceVector,potential.cutoffDistance());
        potential.twoParticleAction(&atom1, &atom2);
        CHECK_ARRAY_CLOSE(forceVector,atom1.force,3,1e-10);
    }

    TEST(USCSiO2) {
        Atom silicon1(AtomType::atomTypeFromAtomType(AtomTypes::Silicon));
        Atom silicon2(AtomType::atomTypeFromAtomType(AtomTypes::Silicon));
        Atom oxygen1(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
        Atom oxygen2(AtomType::atomTypeFromAtomType(AtomTypes::Oxygen));
        vector<Atom*> atoms;
        atoms.push_back(&silicon1);
        atoms.push_back(&silicon2);
        atoms.push_back(&oxygen1);
        atoms.push_back(&oxygen2);

        USCSIO2Potential potential;
        double rMin = 2;
        double rMax = 6;
        int nPoints = 1000;
        double deltaR = (rMax - rMin)/(nPoints - 1);

        vector<double> distances(nPoints);
        vector<double> SiSiForces(nPoints);
        vector<double> SiOForces(nPoints);
        vector<double> OOForces(nPoints);
        for(int i=0; i<nPoints; i++) {
            double x = rMin + i*deltaR;
            distances.at(i) = x;

            // O-O
            for(Atom *atom : atoms) atom->resetForce();
            oxygen2.setPosition(UnitConverter::lengthFromAngstroms(x), 0, 0);
            potential.twoParticleAction(&oxygen1, &oxygen2);
            OOForces.at(i) = oxygen1.force[0];

            // Si-Si
            for(Atom *atom : atoms) atom->resetForce();
            silicon2.setPosition(UnitConverter::lengthFromAngstroms(x), 0, 0);
            potential.twoParticleAction(&silicon1, &silicon2);
            SiSiForces.at(i) = silicon1.force[0];

            // Si-O
            for(Atom *atom : atoms) atom->resetForce();
            oxygen2.setPosition(UnitConverter::lengthFromAngstroms(x), 0, 0);
            potential.twoParticleAction(&silicon1, &oxygen2);
            SiOForces.at(i) = silicon1.force[0];
        }

        cout << "r_MD = [";
        for(int i=0; i<distances.size(); i++) {
            double r = distances.at(i);
            if(i==distances.size()-1) cout << r << "];" << endl;
            else cout << r << ", ";
        }

        cout << "F_Si_Si_MD = [";
        for(int i=0; i<SiSiForces.size(); i++) {
            double F = SiSiForces.at(i);
            if(i==SiSiForces.size()-1) cout << F << "];" << endl;
            else cout << F << ", ";
        }

        cout << "F_Si_O_MD = [";
        for(int i=0; i<SiOForces.size(); i++) {
            double F = SiOForces.at(i);
            if(i==SiOForces.size()-1) cout << F << "];" << endl;
            else cout << F << ", ";
        }

        cout << "F_O_O_MD = [";
        for(int i=0; i<OOForces.size(); i++) {
            double F = OOForces.at(i);
            if(i==OOForces.size()-1) cout << F << "];" << endl;
            else cout << F << ", ";
        }
    }


}
