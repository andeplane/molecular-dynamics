#ifndef ATOM_H
#define ATOM_H
#include <string>
#include <iostream>
#include <atomtype.h>
#include <functional>
#include <random>

using std::string; using std::function;

class AtomType;

class Atom
{
private:
    friend std::ostream& operator<<(std::ostream&stream, const Atom&atom);
    AtomType *m_type;
    unsigned long  m_id;
    unsigned long  m_uniqueId;
    unsigned long  m_originalUniqueId;
    static unsigned long numberOfCreatedAtoms;
    bool m_removed;
    bool m_ghost;
    vector<function<void()>> m_onRemoved;
    vector<function<void()>> m_onMoved;
public:
    double position[3];
    double initial_position[3];
    double velocity[3];
    double force[3];

    Atom();
    Atom(AtomType *atomType);

    void resetForce();
    void move(const double &timestep);
    void kick(const double &timestep, const double oneOverMass);
    void resetVelocityMaxwellian(double temperature);

    AtomType *type() const;
    void setType(AtomType *type);
    bool removed() const;
    void setRemoved(bool removed);
    unsigned long  id() const;
    void setId(unsigned long  id);
    bool ghost() const;
    void setGhost(bool ghost);
    inline void setPosition(const double x, const double y, const double z) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        for(function<void()> onMoved : m_onMoved) {
            onMoved();
        }
    }

    inline void setVelocity(const double x, const double y, const double z) {
        velocity[0] = x;
        velocity[1] = y;
        velocity[2] = z;
    }

    void addOnMoved(const function<void ()> &value);
    void addOnRemoved(const function<void ()> &value);
    unsigned long  uniqueId() const;
    unsigned long originalUniqueId() const;
    void setOriginalUniqueId(unsigned long originalUniqueId);
};

#endif // ATOM_H
