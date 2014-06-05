#include <node.h>
#include <iostream>

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
        if(child->stepIndex() == m_stepIndex) child->step();
    }

    action();
    m_stepIndex++;
}

int Node::frequency() const
{
    return m_frequency;
}

void Node::setFrequency(int frequency)
{
    m_frequency = frequency;
}

Node::Node() :
    m_stepIndex(0),
    m_frequency(1)
{

}
