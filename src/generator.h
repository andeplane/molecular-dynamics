#ifndef GENERATOR_H
#define GENERATOR_H
class System;

class Generator
{
public:
    Generator();
    void generateFCC(System &system, double latticeConstant, int numberOfUnitCells[3], int atomType = 1);
};

#endif // GENERATOR_H
