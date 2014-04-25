#include <iostream>
#include <simulator.h>

using namespace std;

int main()
{
    int numNodesVector[3];
    double systemLength[3];
    numNodesVector[0] = 1; numNodesVector[1] = 1; numNodesVector[2] = 1;
    systemLength[0] = 1; systemLength[1] = 1; systemLength[2] = 1;
    Simulator simulator(0,numNodesVector, systemLength, 2.5);

    return 0;
}

