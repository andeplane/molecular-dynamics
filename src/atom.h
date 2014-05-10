#pragma once
#include <string>
#include <iostream>
#include <atomtype.h>
#include <functional>
#include <random>

using std::string; using std::function;
typedef unsigned long atomUniqueId;

class AtomType;

class Atom
{
private:
    friend std::ostream& operator<<(std::ostream&stream, const Atom&atom);
    AtomType *m_type;
    unsigned long  m_id;
    atomUniqueId  m_uniqueId;
    atomUniqueId  m_originalUniqueId;
    static unsigned long numberOfCreatedAtoms;
    bool m_removed;
    bool m_ghost;
    vector<function<void()>> m_onRemoved;
    vector<function<void()>> m_onMoved;
    vector<unsigned long> m_neighbors;
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

    inline void setPosition(const vector<double> &pos) {
        setPosition(pos.at(0), pos.at(1), pos.at(2));
    }

    inline void setPosition(const double x, const double y, const double z) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        for(function<void()> onMoved : m_onMoved) {
            onMoved();
        }
        m_neighbors.clear(); // This is not valid anymore
    }

    inline void setVelocity(vector<double> vel) {
        setVelocity(vel.at(0), vel.at(1), vel.at(2));
    }

    inline void setVelocity(const double x, const double y, const double z) {
        velocity[0] = x;
        velocity[1] = y;
        velocity[2] = z;
    }

    inline void addVelocity(vector<double> vel) {
        addVelocity(vel.at(0), vel.at(1), vel.at(2));
    }

    inline void addVelocity(const double x, const double y, const double z) {
        setVelocity(velocity[0]+x, velocity[1]+y, velocity[2]+z);
    }

    void addOnMoved(const function<void ()> &value);
    void addOnRemoved(const function<void ()> &value);
    atomUniqueId  uniqueId() const;
    atomUniqueId originalUniqueId() const;
    void setOriginalUniqueId(atomUniqueId originalUniqueId);
    vector<atomUniqueId> &neighbors();
    void addNeighbor(Atom &atom);
};
