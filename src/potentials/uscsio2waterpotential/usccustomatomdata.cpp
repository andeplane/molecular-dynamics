#include <potentials/uscsio2waterpotential/usccustomatomdata.h>
#include <unitconverter.h>
#include <atom.h>

double theta(double r, bool hydrogen) {
    double R = hydrogen ? UnitConverter::lengthFromAngstroms(0.3) : UnitConverter::lengthFromAngstroms(2.0);
    double D = hydrogen ? UnitConverter::lengthFromAngstroms(0.3) : UnitConverter::lengthFromAngstroms(1.4);
    double DMinusR = D-R;
    double DPlusR = D+R;
    if(r <= DMinusR) return 1;
    if(DPlusR < r) return 0;
    else return 1 - (r - R + D)/(2*D) + sin(M_PI*(r - R + D)/D)/(2*M_PI);
}

USCCustomAtomData::CustomAtomData() {

}

void USCCustomAtomData::reset() {
    m_isCalculated = false;
}

void USCCustomAtomData::calculate()
{
    m_nHydrogen = 0;
    m_nSilicon = 0;
    for(Atom *atom : m_atom->neighbors()) {
        double dx = m_atom->position[0] - atom->position[0];
        double dy = m_atom->position[1] - atom->position[1];
        double dz = m_atom->position[2] - atom->position[2];

    }
}
