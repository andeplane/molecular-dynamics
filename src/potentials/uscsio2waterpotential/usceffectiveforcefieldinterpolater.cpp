#include <potentials/uscsio2waterpotential/usceffectiveforcefieldinterpolater.h>
#include <particles/atom.h>
#include <cmath>

void USCEffectiveForceFieldInterpolater::setNSi(double nSi)
{
    m_nSi = nSi;
}

void USCEffectiveForceFieldInterpolater::setNH(double nH)
{
    m_nH = nH;
}
USCEffectiveForceFieldInterpolater::USCEffectiveForceFieldInterpolater() :
    m_nSi(0),
    m_nH(0)
{

}

inline void USCEffectiveForceFieldInterpolater::theta(double r, double R, double D) {
    if(r < R - D) return 1;
    if(r > R + D) return 0;
    else return 1 - (r - R + D)/(2*D) + sin(M_PI*(r - R + D)/D)/(2*M_PI);
}

void USCEffectiveForceFieldInterpolater::compute(Atom *atom)
{
    m_nSi = 0;
    m_nH = 0;
    for(Atom *neighbor : atom->neighbors()) {
        double dr = (atom->position - neighbor->position).length();

        if(neighbor->type().get()->atomicNumber() == AtomTypes::Silicon) {
            m_nSi += theta(dr, m_RSi, m_DSi);
        } else if(neighbor->type().get()->atomicNumber() == AtomTypes::Hydrogen) {
            m_nH += theta(dr, m_RH, m_DH);
        }
    }
}
