#include <cell.h>
#include <atom.h>

Cell::Cell()
{

}

void Cell::addAtom(Atom *atom)
{
    m_atoms.push_back(atom);
}

void Cell::reset() {
    m_atoms.clear();
}

vector<Atom *> Cell::atoms() const
{
    return m_atoms;
}

void Cell::setAtoms(const vector<Atom *> &atoms)
{
    m_atoms = atoms;
}
