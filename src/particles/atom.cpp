#include <particles/atom.h>
#include <string>
#include <random.h>
#include <unitconverter.h>
#include <cmath>
unsigned long  Atom::numberOfCreatedAtoms = 0;

void Atom::addOnRemoved(const function<void ()> &value)
{
    m_onRemoved.push_back(value);
}

void Atom::addOnMoved(const function<void ()> &value)
{
    m_onMoved.push_back(value);
}

void Atom::resetVelocityMaxwellian(double temperature)
{
    if(m_type->atomicNumber() == 0) {
        std::cout << "Warning: Tried to reset maxwellian velocity on an atom of type NoAtom." << std::endl;
    }
    double standardDeviation = sqrt(temperature*m_type->massInverse());
    velocity[0] = Random::nextGauss(0, standardDeviation);
    velocity[1] = Random::nextGauss(0, standardDeviation);
    velocity[2] = Random::nextGauss(0, standardDeviation);
}

void Atom::setOriginalUniqueId(atomUniqueId originalUniqueId)
{
    m_originalUniqueId = originalUniqueId;
}

Atom::Atom() :
    Atom(AtomType::atomTypeFromAtomType(AtomTypes::NoAtom))
{
    
}

Atom::Atom(shared_ptr<AtomType> type) :
    m_type(type),
    m_id(Atom::numberOfCreatedAtoms++),
    m_uniqueId(m_id),
    m_originalUniqueId(m_uniqueId),
    m_removed(false),
    m_ghost(false),
    m_onRemoved(0),
    m_onMoved(0)
{
    memset(position,0,3*sizeof(double));
    memset(velocity,0,3*sizeof(double));
    memset(force,0,3*sizeof(double));
}

unsigned long  Atom::id() const
{
    return m_id;
}

void Atom::setId(unsigned long id)
{
    m_id = id;
}

void Atom::setGhost(bool ghost)
{
    m_ghost = ghost;
}

void Atom::setType(shared_ptr<AtomType> type)
{
    m_type = type;
}

void Atom::move(const double &timestep)
{
    double x = position[0] + velocity[0]*timestep;
    double y = position[1] + velocity[1]*timestep;
    double z = position[2] + velocity[2]*timestep;
    setPosition(x,y,z);
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
        for(function<void()> onRemoved : m_onRemoved) {
            onRemoved();
        }
        m_removed = removed;
    }
}

void Atom::resetForce() {
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
}


std::ostream& operator<<(std::ostream &stream, const Atom &atom) {
    if(!atom.m_ghost) return stream << "Atom (unique id " << atom.uniqueId() << ", original unique id " << atom.originalUniqueId() << ") of type '" << atom.type()->name() << "' with id " << atom.id() << " r" << atom.uniqueId() << "=[" << UnitConverter::lengthToAngstroms(atom.position[0]) << ", " << UnitConverter::lengthToAngstroms(atom.position[1]) << ", " << UnitConverter::lengthToAngstroms(atom.position[2]) << "]  v" << atom.uniqueId() << "=[" << atom.velocity[0] << ", " << atom.velocity[1] << ", " << atom.velocity[2] << "]  f" << atom.uniqueId() << "=[" << atom.force[0] << ", " << atom.force[1] << ", " << atom.force[2] << "]";
    else return stream << "Ghost atom (unique id " << atom.uniqueId() << ", original unique id " << atom.originalUniqueId() << ") of type '" << atom.type()->name() << "' with id " << atom.id() << " r" << atom.uniqueId() << "=[" << UnitConverter::lengthToAngstroms(atom.position[0]) << ", " << UnitConverter::lengthToAngstroms(atom.position[1]) << ", " << UnitConverter::lengthToAngstroms(atom.position[2]) << "]  v" << atom.uniqueId() << "=[" << atom.velocity[0] << ", " << atom.velocity[1] << ", " << atom.velocity[2] << "]  f" << atom.uniqueId() << "=[" << atom.force[0] << ", " << atom.force[1] << ", " << atom.force[2] << "]";

}
