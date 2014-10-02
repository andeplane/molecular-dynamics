#include <cell.h>
#include <particles/atom.h>


int Cell::cellIndex() const
{
    return m_cellIndex;
}

void Cell::setCellIndex(int cellIndex)
{
    m_cellIndex = cellIndex;
}
Cell::Cell() :
    m_cellIndex(-1)
{

}

void Cell::addAtom(Atom *atom)
{
    m_atoms.push_back(atom);
}

void Cell::reset() {
    m_atoms.clear();
}

vector<Atom *> &Cell::atoms()
{
    return m_atoms;
}

void Cell::setAtoms(const vector<Atom *> &atoms)
{
    m_atoms = atoms;
}
