#pragma once
#include <memory>
#include <node.h>
#include <atom.h>
#include <unordered_map>
using std::unordered_map;
using std::shared_ptr;

class System;

class NeighborList : public Node
{
private:
    shared_ptr<System> m_system;
    unordered_map<atomUniqueId, shared_ptr<Atom>> m_list;
    double m_maxDistance;

public:
    NeighborList(shared_ptr<System> system, double maxDistance);
    double maxDistance() const;
    void setMaxDistance(double maxDistance);
    virtual void action();
    unordered_map<atomUniqueId, shared_ptr<Atom> > list() const;
    void setList(const unordered_map<atomUniqueId, shared_ptr<Atom> > &list);
};
