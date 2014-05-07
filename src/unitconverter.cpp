#include <unitconverter.h>
#include <iostream>
#include <cmath>
using namespace std;

UnitConverter::UnitConverter()
{
    // Atomic units
    // [1] http://en.wikipedia.org/wiki/Atomic_units

    m0 = 9.10938291e-31;  // SI [kg]
    q0 = 1.602176565e-19; // SI [C]
    hbar0 = 1.054571726e-34; // SI [Js]
    electricConstant0 = 8.9875517873681e9; // SI [kgm^3/(s^-2 C^-2)]
    kb = 1.3806488e-23; // SI [J/K]

    a0 = hbar0*hbar0/(electricConstant0*m0*q0*q0);
    E0 = m0*q0*q0*q0*q0*electricConstant0*electricConstant0/(hbar0*hbar0);
    t0 = hbar0/E0;
    v0 = a0*E0/hbar0;
    F0 = E0/a0;
    T0 = E0/kb;
    P0 = E0/(a0*a0*a0);
    visc0 = P0*t0;
    diff0 = a0*a0/t0;
}

double UnitConverter::pressureToSI(double P) { return P0*P; }
double UnitConverter::pressureFromSI(double P) { return P/P0; }

double UnitConverter::temperatureToSI(double T) { return T0*T; }
double UnitConverter::temperatureFromSI(double T) { return T/T0; }

double UnitConverter::massToSI(double m) { return m0*m; }
double UnitConverter::massFromSI(double m) { return m/m0; }

double UnitConverter::chargeToSI(double q) { return q*q0; }
double UnitConverter::chargeFromSI(double q) { return q/q0; }

double UnitConverter::hbarToSI(double hbar) { return hbar*hbar0; }
double UnitConverter::hbarFromSI(double hbar) { return hbar/hbar0; }

double UnitConverter::electricConstantToSI(double electricConstant) { return electricConstant*electricConstant0; }
double UnitConverter::electricConstantFromSI(double electricConstant) { return electricConstant/electricConstant0; }

double UnitConverter::lengthToSI(double L) { return a0*L; }
double UnitConverter::lengthFromSI(double L) { return L/a0; }

double UnitConverter::lengthToAngstroms(double L) { return a0*L*1e10; }
double UnitConverter::lengthFromAngstroms(double L) { return L/(a0*1e10); }

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


