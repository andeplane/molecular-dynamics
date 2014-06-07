#include <statistics/paircorrelationsampler.h>
#include <system.h>
#include <statistics/neighborlist.h>

int PairCorrelationSampler::numberOfBins() const
{
    return m_numberOfBins;
}

void PairCorrelationSampler::setNumberOfBins(int numberOfBins)
{
    m_numberOfBins = numberOfBins;
}
PairCorrelationSampler::PairCorrelationSampler(shared_ptr<AtomType> atomType1, shared_ptr<AtomType> atomType2, shared_ptr<System> system, shared_ptr<NeighborList> neighborList, int numberOfBins) :
    m_atomType1(atomType1),
    m_atomType2(atomType2),
    m_system(system),
    m_numberOfBins(numberOfBins),
    m_neighborList(neighborList)
{
    addInput(neighborList);
}

void PairCorrelationSampler::action()
{
    vector<int> pairCorrelationFunction(m_numberOfBins, 0);

    m_system.lock()->atomManager().atoms().iterate([&] (Atom &atom) {
        if(atom.type() != m_atomType1) return;

        atomUniqueId uniqueId = atom.uniqueId();
        auto &neighborsOfThisAtom = m_neighborList.lock()->list()[uniqueId];
        for(auto neighborAtom : neighborsOfThisAtom) {
            if(neighborAtom->type() != m_atomType2) continue;
            double dx = atom.position[0] - neighborAtom->position[0];
            double dy = atom.position[1] - neighborAtom->position[1];
            double dz = atom.position[2] - neighborAtom->position[2];
            double dr = sqrt(dx*dx + dy*dy + dz*dz);
            int binIndex = dr / maxDistance()*m_numberOfBins;
            if(binIndex < m_numberOfBins) {
                pairCorrelationFunction.at(binIndex)++;
            }
        }
    });

    m_value.addValue(pairCorrelationFunction);
    pairCorrelationFunction.clear();
}

double PairCorrelationSampler::maxDistance()
{
    return m_neighborList.lock()->maxDistance();
}
