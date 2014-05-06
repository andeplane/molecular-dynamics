#pragma once

class UnitConverter
{
public:
    double m0;
    double L0;
    double L0_angstrom;
    double E0;
    double E0ev;
    double e;
    double kb;
    double t0;
    double F0;
    double T0;
    double P0;
    double v0;
    double visc0;
    double diff0;

    UnitConverter();

    double pressureToSI(double P);
    double pressureFromSI(double P);

    double temperatureToSI(double T);
    double temperatureFromSI(double T);

    double massToSI(double m);
    double massFromSI(double m);

    double lengthToSI(double L);
    double lengthFromSI(double L);

    double lengthToAngstroms(double L);
    double lengthFromAngstroms(double L);

    double forceToSI(double F);
    double forceFromSI(double F);

    double energyToSI(double E);
    double energyFromSI(double E);

    double timeToSI(double t);
    double timeFromSI(double t);

    double velocityToSI(double v);
    double velocityFromSI(double v);

    double viscosityToSI(double v);
    double viscosityFromSI(double v);

    double diffusionToSI(double d);
    double diffusionFromSI(double d);

    double energyToEv(double E);
    double energyFromEv(double E);

    double degreesToRadians(double v);
    double radiansToDegrees(double v);
    double chargeToSI(double m);
    double chargeFromSI(double m);
};
