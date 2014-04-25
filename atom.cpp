#include "atom.h"
#include <string>

int Atom::id() const
{
    return m_id;
}

void Atom::setId(int id)
{
    m_id = id;
}
// Default atom is Hydrogen - atom type 1
Atom::Atom() :
    m_type(1),
    m_moved(false),
    m_id(-1)
{
    memset(position,0,3*sizeof(double));
    memset(velocity,0,3*sizeof(double));
    memset(force,0,3*sizeof(double));
}

int Atom::type() const
{
    return m_type;
}

void Atom::setType(int type)
{
    m_type = type;
}

void Atom::move(double &timestep)
{
    position[0] += velocity[0]*timestep;
    position[1] += velocity[1]*timestep;
    position[2] += velocity[2]*timestep;
}

void Atom::kick(double &timestep, double oneOverMass)
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
