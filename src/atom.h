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
    int m_id;
    static int numberOfCreatedAtoms;
    bool m_removed;
    bool m_ghost;
    function<void()> m_onRemoved;
    function<void()> m_onMoved;
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
    void resetMaxwellianVelocity(double temperature);

    AtomType *type() const;
    void setType(AtomType *type);
    bool removed() const;
    void setRemoved(bool removed);
    int id() const;
    void setId(int id);
    bool ghost() const;
    void setGhost(bool ghost);
    void setOnRemoved(const function<void ()> &value);
    void setOnMoved(const function<void ()> &onMoved);
    inline void setPosition(const double x, const double y, const double z) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        if(m_onMoved) m_onMoved();
    }

    inline void setVelocity(const double x, const double y, const double z) {
        velocity[0] = x;
        velocity[1] = y;
        velocity[2] = z;
    }
};

#endif // ATOM_H
