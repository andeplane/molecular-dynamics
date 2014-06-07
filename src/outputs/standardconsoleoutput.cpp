#include <outputs/standardconsoleoutput.h>
#include <unitconverter.h>

int StandardConsoleOutput::frequency() const
{
    return m_frequency;
}

void StandardConsoleOutput::setFrequency(int frequency)
{
    m_frequency = frequency;
}

shared_ptr<StandardConsoleOutput> StandardConsoleOutput::create(shared_ptr<Measurements::TotalEnergySampler> totalEnergySampler, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system, int frequency)
{
    return std::make_shared<StandardConsoleOutput>(totalEnergySampler, temperatureSampler, system, frequency);
}

StandardConsoleOutput::StandardConsoleOutput(shared_ptr<Measurements::TotalEnergySampler> totalEnergySampler, shared_ptr<Measurements::TemperatureSampler> temperatureSampler, shared_ptr<System> system, int frequency) :
    m_totalEnergySampler(totalEnergySampler),
    m_temperatureSampler(temperatureSampler),
    m_system(system),
    m_frequency(frequency)
{
    addInput(m_totalEnergySampler);
    addInput(m_temperatureSampler);
}

void StandardConsoleOutput::action()
{
    if(updatedAtStep() % m_frequency == 0) {
        double energyEv = UnitConverter::energyToEv(m_totalEnergySampler->value().currentValueScalar());
        double energyEvPerAtom = energyEv / std::max(1,m_system->numberOfAtoms());
        double temperatureSI = UnitConverter::temperatureToSI(m_temperatureSampler->value().currentValueScalar());
        cout << updatedAtStep() << ": E/N=" << energyEvPerAtom << " eV   T=" << temperatureSI << " K " << endl;
    }
}
