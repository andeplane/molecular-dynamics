#include "atom.h"
#include <string>

Atom::Atom() :
    m_type(AtomType::atomTypeFromAtomType(AtomTypes::NoAtom)),
    m_moved(false),
    m_id(-1),
    m_ghost(false)
{
    memset(position,0,3*sizeof(double));
    memset(velocity,0,3*sizeof(double));
    memset(force,0,3*sizeof(double));
}

Atom::Atom(AtomType *type) :
    m_type(type)
{
    Atom::Atom(); // Call default constructor
}

int Atom::id() const
{
    return m_id;
}

void Atom::setId(int id)
{
    m_id = id;
}

bool Atom::ghost() const
{
    return m_ghost;
}

void Atom::setGhost(bool ghost)
{
    m_ghost = ghost;
}

AtomType *Atom::type() const
{
    return m_type;
}

void Atom::setType(AtomType *type)
{
    m_type = type;
}

void Atom::move(const double &timestep)
{
    position[0] += velocity[0]*timestep;
    position[1] += velocity[1]*timestep;
    position[2] += velocity[2]*timestep;
}

void Atom::kick(const double &timestep, const double oneOverMass)
{
    velocity[0] += force[0]*timestep*oneOverMass;
    velocity[1] += force[1]*timestep*oneOverMass;
    velocity[2] += force[2]*timestep*oneOverMass;
}


bool Atom::moved() const
{
    return m_moved;
}

void Atom::setMoved(bool moved)
{
    m_moved = moved;
}

void Atom::resetForce() {
    memset(force,0,3*sizeof(double));
}

std::ostream& operator<<(std::ostream &stream, const Atom &atom) {
    return stream << "Atom of type " << atom.type()->name() << " r=(" << atom.position[0] << ", " << atom.position[1] << ", " << atom.position[2] << ")  v=(" << atom.velocity[0] << ", " << atom.velocity[1] << ", " << atom.velocity[2] << ")  f= (" << atom.force[0] << ", " << atom.force[1] << ", " << atom.force[2] << ")";
}
