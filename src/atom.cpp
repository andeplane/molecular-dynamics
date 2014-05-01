#include <atom.h>
#include <string>
#include <random.h>

int Atom::numberOfCreatedAtoms = 0;

void Atom::setOnRemoved(const function<void ()> &value)
{
    m_onRemoved = value;
}

void Atom::setOnMoved(const function<void ()> &onMoved)
{
    m_onMoved = onMoved;
}

void Atom::resetMaxwellianVelocity(double temperature)
{
    double standardDeviation = sqrt(temperature*m_type->massInverse());
    velocity[0] = Random::nextGauss(0, standardDeviation);
    velocity[1] = Random::nextGauss(0, standardDeviation);
    velocity[2] = Random::nextGauss(0, standardDeviation);
}

Atom::Atom() :
    Atom(AtomType::atomTypeFromAtomType(AtomTypes::NoAtom))
{

}

Atom::Atom(AtomType *type) :
    m_id(Atom::numberOfCreatedAtoms++),
    m_removed(false),
    m_ghost(false),
    m_onRemoved(0),
    m_onMoved(0),
    m_type(type)
{
    memset(position,0,3*sizeof(double));
    memset(velocity,0,3*sizeof(double));
    memset(force,0,3*sizeof(double));
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

    m_onMoved();
}

void Atom::kick(const double &timestep, const double oneOverMass)
{
    velocity[0] += force[0]*timestep*oneOverMass;
    velocity[1] += force[1]*timestep*oneOverMass;
    velocity[2] += force[2]*timestep*oneOverMass;
}


bool Atom::removed() const
{
    return m_removed;
}

void Atom::setRemoved(bool removed)
{
    if(m_removed != removed) {
        m_onRemoved();
        m_removed = removed;
    }
}

void Atom::resetForce() {
    memset(force,0,3*sizeof(double));
}

std::ostream& operator<<(std::ostream &stream, const Atom &atom) {
    if(!atom.m_ghost) return stream << "Atom of type '" << atom.type()->name() << "' with id " << atom.id() << " r=(" << atom.position[0] << ", " << atom.position[1] << ", " << atom.position[2] << ")  v=(" << atom.velocity[0] << ", " << atom.velocity[1] << ", " << atom.velocity[2] << ")  f= (" << atom.force[0] << ", " << atom.force[1] << ", " << atom.force[2] << ")";
    else return stream << "Ghost atom of type '" << atom.type()->name() << "' with id " << atom.id() << " r=(" << atom.position[0] << ", " << atom.position[1] << ", " << atom.position[2] << ")  v=(" << atom.velocity[0] << ", " << atom.velocity[1] << ", " << atom.velocity[2] << ")  f= (" << atom.force[0] << ", " << atom.force[1] << ", " << atom.force[2] << ")";

}
