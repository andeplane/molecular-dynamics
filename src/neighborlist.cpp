#include <neighborlist.h>
#include <atom.h>

NeighborMap &NeighborList::neighborMap()
{
    return m_neighborMap;
}

void NeighborList::addNeighbors(Atom *atom1, Atom *atom2)
{
    m_neighborMap[atom1->uniqueId()].push_back(atom2);
    m_neighborMap[atom1->uniqueId()].push_back(atom1);
}

vector<Atom *> &NeighborList::neighborsForAtom(Atom *atom)
{
    return m_neighborMap[atom->uniqueId()];
}

NeighborList::NeighborList()
{

}

NeighborList::~NeighborList()
{
    for(NeighborMap::iterator it = m_neighborMap.begin(); it!= m_neighborMap.end(); ++it) {
        vector<Atom*> &list = (*it).second;
        list.clear();
    }

    m_neighborMap.clear();
}

