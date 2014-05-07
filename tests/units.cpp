#include <vector>
#include <iostream>
#include <cmath>
using std::vector;
using std::cout; using std::endl;

#include <UnitTest++.h>
#include <includes.h>
#include <unitconverter.h>

SUITE(Units) {
    TEST(AtomicUnits) {
        // Units from
        // [1] http://en.wikipedia.org/wiki/Atomic_units

        UnitConverter uc;

        // Fundamental units
        CHECK_CLOSE(9.10938291e-31, uc.massToSI(1.0), 1e-38); // Mass
        CHECK_CLOSE(1.602176565e-19, uc.chargeToSI(1.0), 1e-26); // Charge
        CHECK_CLOSE(1.054571726e-34, uc.hbarToSI(1.0), 1e-40); // hbar
        CHECK_CLOSE(8.9875517873681e9, uc.electricConstantToSI(1.0), 1e2); // Electric constant
        // Derived units
        CHECK_CLOSE(5.2917721092e-11, uc.lengthToSI(1.0), 1e-19); // Length
        CHECK_CLOSE(4.35974417e-18, uc.energyToSI(1.0), 1e-24); // Energy
        CHECK_CLOSE(2.418884326505e-17, uc.timeToSI(1.0), 1e-25); // Time
        CHECK_CLOSE(2.1876912633e6, uc.velocityToSI(1.0), 1e-3); // Velocity
        CHECK_CLOSE(8.2387225e-8, uc.forceToSI(1.0), 1e-14); // Force
        CHECK_CLOSE(3.1577464e5, uc.temperatureToSI(1.0), 1); // Temperature
        CHECK_CLOSE(2.9421912e13, uc.pressureToSI(1.0), 1e9); // Pressure

    }
}
