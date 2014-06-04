#pragma once
#include <vector>
#include <memory>
using std::vector;
using std::shared_ptr;

class Node
{
private:
    vector<shared_ptr<Node> > m_children;
public:
    Node();
    vector<shared_ptr<Node> > children();
    void addChild(shared_ptr<Node> object);
    shared_ptr<Node> getChildByIndex(decltype(m_children.size()) index);
    virtual void action() { }
    void step();
};
