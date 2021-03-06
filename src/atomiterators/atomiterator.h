#pragma once
#include <functional>
#include <vector>
using std::vector;
using std::function;

#include <atommanager.h>
#include <particles/atom.h>

class AtomIterator
{
protected:
    function<void(Atom *atom1, Atom *atom2)> m_twoParticleAction;
    function<void(Atom *atom1, Atom *atom2, Atom *atom3)> m_threeParticleAction;
public:
    AtomIterator();
    virtual void iterate(AtomManager &atomManager) = 0;

    function<void (Atom *atom1, Atom *atom2)> twoParticleAction() const;
    void setTwoParticleAction(const function<void (Atom *atom1, Atom *atom2)> &twoParticleAction);
    function<void (Atom *atom1, Atom *atom2, Atom *atom3)> threeParticleAction() const;
    virtual void setThreeParticleAction(const function<void (Atom *atom1, Atom *atom2, Atom *atom3)> &threeParticleAction);
};
