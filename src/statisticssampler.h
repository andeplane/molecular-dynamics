#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
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

#endif // STATISTICSSAMPLER_H
