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
    velocity.randomGaussian(0, standardDeviation);
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
    position.addAndMultiply(velocity, timestep);
}

void Atom::kick(const double &timestep, const double oneOverMass)
{
    velocity.addAndMultiply(force, timestep*oneOverMass);
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
    force.setToZero();
}


std::ostream& operator<<(std::ostream &stream, const Atom &atom) {
    if(!atom.ghost()) return stream << "Atom (unique id " << atom.uniqueId() << ", original unique id " << atom.originalUniqueId() << ") of type '" << atom.type()->name() << "' with id " << atom.id() << " r" << atom.uniqueId() << "=[" << UnitConverter::lengthToAngstroms(atom.position.x()) << ", " << UnitConverter::lengthToAngstroms(atom.position.y()) << ", " << UnitConverter::lengthToAngstroms(atom.position.z()) << "]  v" << atom.uniqueId() << "=[" << atom.velocity.x() << ", " << atom.velocity.y() << ", " << atom.velocity.z() << "]  f" << atom.uniqueId() << "=[" << atom.force.x() << ", " << atom.force.y() << ", " << atom.force.z() << "]";
    else return stream << "Ghost atom (unique id " << atom.uniqueId() << ", original unique id " << atom.originalUniqueId() << ") of type '" << atom.type()->name() << "' with id " << atom.id() << " r" << atom.uniqueId() << "=[" << UnitConverter::lengthToAngstroms(atom.position.x()) << ", " << UnitConverter::lengthToAngstroms(atom.position.y()) << ", " << UnitConverter::lengthToAngstroms(atom.position.z()) << "]  v" << atom.uniqueId() << "=[" << atom.velocity.x() << ", " << atom.velocity.y() << ", " << atom.velocity.z() << "]  f" << atom.uniqueId() << "=[" << atom.force.x() << ", " << atom.force.y() << ", " << atom.force.z() << "]";

}
