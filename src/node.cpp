#include <node.h>
void Node::step(int stepIndex)
{
    for(auto input : m_inputs) { input->step(stepIndex); }
    if(stepIndex != m_updatedAtStep) {
        m_updatedAtStep = stepIndex;
        action();
    }
}

vector<shared_ptr<Node> > Node::inputs() const
{
    return m_inputs;
}

void Node::addInput(shared_ptr<Node> input)
{
    m_inputs.push_back(input);
}


int Node::updatedAtStep() const
{
    return m_updatedAtStep;
}

Node::Node() :
    m_updatedAtStep(-1)
{

}
