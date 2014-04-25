#include "systemcell.h"

#include <atom.h>

int SystemCell::numberOfCells[3];
double SystemCell::cellLength[3];

SystemCell::SystemCell()
{

}

void SystemCell::addAtom(Atom *atom)
{
    m_atoms.push_back(atom);
}

void SystemCell::reset() {
    m_atoms.clear();
}

vector<Atom *> SystemCell::atoms() const
{
    return m_atoms;
}

void SystemCell::setAtoms(const vector<Atom *> &atoms)
{
    m_atoms = atoms;
}
