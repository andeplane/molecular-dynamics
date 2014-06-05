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

set<shared_ptr<Node>> Node::getChildByTag(int tag, bool recursive)
{
    set<shared_ptr<Node>> children;

    for(auto child : m_children) {
        if(child->tag() == tag) {
            children.insert(child);
            if(recursive) {
                auto childrenOfChild = child->getChildByTag(tag, recursive);
                children.insert(childrenOfChild.begin(), childrenOfChild.end());
            }
        }
    }

    return children;
}

void Node::step()
{
    for(auto child : m_children) {
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

int Node::tag() const {
    return m_tag;
}

void Node::setTag(int tag) {
    m_tag = tag;
}

Node::Node() :
    m_stepIndex(0),
    m_frequency(1),
    m_tag(0)
{

}
