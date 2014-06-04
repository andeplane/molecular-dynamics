#include <mdnode.h>

vector<shared_ptr<MDNode> > MDNode::children()
{
    return m_children;
}

void MDNode::addChild(shared_ptr<MDNode> object)
{
    m_children.push_back(object);
}

shared_ptr<MDNode> MDNode::getChildByIndex(decltype(m_children.size()) index)
{
    return m_children.at(index);
}

void MDNode::step()
{
    for(shared_ptr<MDNode> &child : m_children) {
        child->step();
    }

    action();
}

MDNode::MDNode()
{

}
