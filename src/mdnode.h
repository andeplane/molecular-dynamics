#pragma once
#include <vector>
#include <memory>
using std::vector;
using std::shared_ptr;

class MDNode
{
private:
    vector<shared_ptr<MDNode> > m_children;
public:
    MDNode();
    vector<shared_ptr<MDNode> > children();
    void addChild(shared_ptr<MDNode> object);
    shared_ptr<MDNode> getChildByIndex(decltype(m_children.size()) index);
    virtual void action() = 0;
    void step();
};
