#pragma once
#include <vector>
using std::vector;

class System;
class StatisticsSampler
{
public:
    StatisticsSampler();
    double calculateKineticEnergy(System &system);
    double calculatePotentialEnergy(System &system);
    double calculateTotalEnergy(System &system);
    vector<double> calculateTotalMomentum(System &system);
};
