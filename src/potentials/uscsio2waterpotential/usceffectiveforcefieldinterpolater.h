#pragma once
#include <particles/customproperty.h>
#include <unitconverter.h>
class Atom;

class USCEffectiveForceFieldInterpolater : public CustomProperty
{
private:
    const double m_RSi = UnitConverter::lengthFromAngstroms(2.0);
    const double m_DSi = UnitConverter::lengthFromAngstroms(0.3);
    const double m_RH = UnitConverter::lengthFromAngstroms(1.4);
    const double m_DH = UnitConverter::lengthFromAngstroms(0.3);

    double m_nSi;
    double m_nH;

    void theta(double r, double R, double D);
public:
    USCEffectiveForceFieldInterpolater();
    void compute(Atom *atom);

    inline double nSi() const { return m_nSi; }
    void setNSi(double nSi);
    inline double nH() const { return m_nH; }
    void setNH(double nH);
};
