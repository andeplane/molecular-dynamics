#include "random.h"

std::random_device Random::rd;
std::mt19937 Random::e2;
bool Random::initialized = false;

void Random::initialize()
{
    Random::initialized = true;
    Random::e2 = std::mt19937(Random::rd());
}

Random::Random()
{
}

double Random::nextGauss(const double &mean, const double &standardDeviation)
{
    std::normal_distribution<> normal_dist(mean, standardDeviation);
    return normal_dist(Random::e2);
}

void Random::setSeed(int seed)
{
    Random::e2.seed(seed);
}
