#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
#include <unordered_map>
#include <vector>
using std::unordered_map; using std::vector;
class Atom;

typedef unordered_map<unsigned long,vector<Atom*> > NeighborMap;
class NeighborList
{
private:
    NeighborMap m_neighborMap;
public:
    NeighborList();
    ~NeighborList();
    NeighborMap &neighborMap();
    void addNeighbors(Atom *atom1, Atom *atom2);
    vector<Atom*> &neighborsForAtom(Atom *atom);
};

#endif // NEIGHBORLIST_H
