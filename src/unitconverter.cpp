#include <unitconverter.h>
#include <iostream>
#include <cmath>
using namespace std;

UnitConverter::UnitConverter()
{
    m0 = 1.66053886e-27;  // SI [kg]
    L0 = 3.405e-10;       // SI [m]
    L0_angstrom = 3.405;  // Angstroms
    E0 = 1.65088e-21;     // SI [J]
    E0ev = 1.030398291e-2;
    kb = 1.3806503e-23;   // SI [J/K]
    e  = 1.60217657e-19;  // SI [C]

    t0 = L0*sqrt(m0/E0);
    F0 = E0/L0;
    T0 = E0/kb;
    P0 = m0/(t0*t0*L0);
    v0 = L0/t0;
    visc0 = P0*t0;
    diff0 = L0*L0/t0;
}

double UnitConverter::pressureToSI(double P) { return P0*P; }
double UnitConverter::pressureFromSI(double P) { return P/P0; }

double UnitConverter::temperatureToSI(double T) { return T0*T; }
double UnitConverter::temperatureFromSI(double T) { return T/T0; }

double UnitConverter::massToSI(double m) { return m0*m; }
double UnitConverter::massFromSI(double m) { return m/m0; }

double UnitConverter::chargeToSI(double q) { return 0; }
double UnitConverter::chargeFromSI(double q) { return 0; }

double UnitConverter::lengthToSI(double L) { return L0*L; }
double UnitConverter::lengthFromSI(double L) { return L/L0; }

double UnitConverter::lengthToAngstroms(double L) { return L0_angstrom*L; }
double UnitConverter::lengthFromAngstroms(double L) { return L/L0_angstrom; }

double UnitConverter::forceToSI(double F) { return F0*F; }
double UnitConverter::forceFromSI(double F) { return F/F0; }

double UnitConverter::energyToSI(double E) { return E0*E; }
double UnitConverter::energyFromSI(double E) { return E/E0; }

double UnitConverter::energyToEv(double E) { return E0*E; }
double UnitConverter::energyFromEv(double E) { return E/E0; }

double UnitConverter::degreesToRadians(double v) { return M_PI/180*v; }
double UnitConverter::radiansToDegrees(double v) { return 180/M_PI*v; }

double UnitConverter::timeToSI(double t) { return t0*t; }
double UnitConverter::timeFromSI(double t) { return t/t0; }

double UnitConverter::velocityToSI(double v) { return v*v0; }
double UnitConverter::velocityFromSI(double v) { return v/v0; }

double UnitConverter::viscosityToSI(double v) { return v*visc0; }
double UnitConverter::viscosityFromSI(double v) { return v/visc0; }

double UnitConverter::diffusionToSI(double d) { return d*diff0; }
double UnitConverter::diffusionFromSI(double d) { return d/diff0; }


