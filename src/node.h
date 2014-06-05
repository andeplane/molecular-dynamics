#pragma once
#include <vector>
#include <memory>
#include <typeinfo>
#include <set>
using std::vector;
using std::set;
using std::shared_ptr;

//template <class T>
//class NodeList {
//private:
//    vector<T> m_nodes
//public:
//    NodeList()
//};

class Node
{
private:
    shared_ptr<Node> m_parent;
    vector<shared_ptr<Node> > m_children;
    int m_stepIndex;
public:
    Node();
    vector<shared_ptr<Node> > children();
    void addChild(shared_ptr<Node> object);
    shared_ptr<Node> getChildByIndex(decltype(m_children.size()) index);
    virtual void action() { }
    void step();
    inline int stepIndex() { return m_stepIndex; }

    template <class T>
    set<shared_ptr<T> > getChildrenOfClass(bool recursive) {
        set<shared_ptr<T> > children;
        for(auto child : m_children) {
            auto castedChild = std::dynamic_pointer_cast<T>(child);
            if(castedChild != 0) {
                children.insert(castedChild);
            }

            if(recursive) {
                auto recursiveChildren = child->getChildrenOfClass<T>(true);
                children.insert(recursiveChildren.begin(), recursiveChildren.end());
            }
        }

        return children;
    }
};
