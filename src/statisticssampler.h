#pragma once
#include <vector>
using std::vector;
using std::shared_ptr;

class System;
class StatisticsSampler
{
public:
    StatisticsSampler();
    double calculateKineticEnergy(shared_ptr<System> system);
    double calculatePotentialEnergy(shared_ptr<System> system);
    double calculateTotalEnergy(shared_ptr<System> system);
    vector<double> calculateTotalMomentum(shared_ptr<System> system);
    double calculateTemperature(shared_ptr<System> system);
};
