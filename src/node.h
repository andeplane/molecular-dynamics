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
    vector<shared_ptr<Node> > m_inputs;
    int m_updatedAtStep;
public:
    Node();
    // Getters/setters
    vector<shared_ptr<Node> > inputs() const;
    void addInput(shared_ptr<Node> input);

    // Action functions
    virtual void action() { }
    void step(int stepIndex);
    int updatedAtStep() const;
};
