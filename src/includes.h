#ifndef INCLUDES_H
#define INCLUDES_H
#include <system.h>
#include <simulator.h>
#include <potentials/potential.h>
#include <integrators/integrator.h>
#include <atom.h>
#include <generator.h>
#include <atomlist.h>
#include <atommanager.h>
#include <potentials/potentials.h>
#include <integrators/integrators.h>
#include <atomiterators/atomiterators.h>
#include <vector>

std::ostream& operator<<(std::ostream &stream, std::vector<double> &vec);

#endif // INCLUDES_H
