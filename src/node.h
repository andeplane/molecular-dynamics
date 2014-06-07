#pragma once
#include <vector>
#include <memory>
#include <typeinfo>
#include <set>
using std::vector;
using std::set;
using std::shared_ptr;
using std::weak_ptr;

class Node
{
private:
    vector<shared_ptr<Node> > m_children;
    vector<weak_ptr<Node> > m_dependencies;
    int m_stepIndex;
    int m_frequency;
    int m_tag;
public:
    Node();
    // Getters/setters
    int frequency() const;
    void setFrequency(int frequency);
    int tag() const;
    void setTag(int tag);
    inline int stepIndex() { return m_stepIndex; }
    vector<weak_ptr<Node> > dependencies() const;
    void addDependency(weak_ptr<Node> dependency);

    // Action functions
    virtual void action() { }
    void step(int stepIndex);
    void step();
    // Children
    vector<shared_ptr<Node> > children();
    void addChild(shared_ptr<Node> object);
    shared_ptr<Node> getChildByIndex(decltype(m_children.size()) index);
    set<shared_ptr<Node>> getChildByTag(int tag, bool recursive = true);
    template <class T>
    set<shared_ptr<T> > getChildrenOfClass(bool recursive = true) {
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
