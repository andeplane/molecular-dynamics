#ifndef ATOM_H
#define ATOM_H

class Atom
{
private:
    int m_id;
    int m_type;
    bool m_moved;
    bool m_isGhost;
public:
    Atom();
    void resetForce();
    int type() const;
    void setType(int type);
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
};

#endif // ATOM_H
