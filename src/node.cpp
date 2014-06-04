#include <node.h>

vector<shared_ptr<Node> > Node::children()
{
    return m_children;
}

void Node::addChild(shared_ptr<Node> object)
{
    m_children.push_back(object);
}

shared_ptr<Node> Node::getChildByIndex(decltype(m_children.size()) index)
{
    return m_children.at(index);
}

void Node::step()
{
    for(shared_ptr<Node> &child : m_children) {
        child->step();
    }

    action();
}

Node::Node()
{

}
