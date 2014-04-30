#ifndef RANDOM_H
#define RANDOM_H
#include <random>

class Random
{
private:
    static std::random_device rd;
    static std::mt19937 e2;
    static bool initialized;
    static void initialize();
public:
    Random();
    static double nextGauss(const double &mean, const double &standardDeviation);
    static void setSeed(int seed);
};

#endif // RANDOM_H
