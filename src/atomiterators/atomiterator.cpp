#include "atomiterator.h"

function<void (Atom *atom1, Atom *atom2)> AtomIterator::twoParticleAction() const
{
    return m_twoParticleAction;
}

void AtomIterator::setTwoParticleAction(const function<void (Atom *atom1, Atom *atom2)> &twoParticleAction)
{
    m_twoParticleAction = twoParticleAction;
}

function<void (Atom *atom1, Atom *atom2, Atom *atom3)> AtomIterator::threeParticleAction() const
{
    return m_threeParticleAction;
}

void AtomIterator::setThreeParticleAction(const function<void (Atom *atom1, Atom *atom2, Atom *atom3)> &threeParticleAction)
{
    m_threeParticleAction = threeParticleAction;
}

AtomIterator::AtomIterator() :
    m_twoParticleAction(0),
    m_threeParticleAction(0)
{

}
