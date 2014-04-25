#include "systemcell.h"

#include <atom.h>

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
