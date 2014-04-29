#ifndef ATOM_H
#define ATOM_H
#include <string>
#include <iostream>
#include <atomtype.h>

using std::string;

class AtomType;

class Atom
{
private:
    friend std::ostream& operator<<(std::ostream&stream, const Atom&atom);
    AtomType *m_type;
    int m_id;
    bool m_moved;
    bool m_isGhost;
public:
    Atom();
    Atom(AtomType *atomType);
    void resetForce();
    void move(double &timestep);
    void kick(double &timestep, double oneOverMass);

    double position[3];
    double initial_position[3];
    double velocity[3];
    double force[3];
    bool moved() const;
    void setMoved(bool moved);
    int id() const;
    void setId(int id);
    bool isGhost() const;
    void setIsGhost(bool isGhost);
    AtomType *type() const;
    void setType(AtomType *type);
};

#endif // ATOM_H
