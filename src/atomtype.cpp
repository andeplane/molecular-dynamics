#include "atomtype.h"

bool AtomType::isInitialized = false;
vector<AtomType*> AtomType::atomTypes;

AtomType::AtomType(int atomicNumber, double mass, string name, string symbol) :
    m_atomicNumber(atomicNumber),
    m_mass(mass),
    m_name(name),
    m_symbol(symbol),
    m_massInverse(1.0/mass)
{

}

AtomType *AtomType::atomTypeFromAtomicNumber(int atomicNumber) {
    if(!AtomType::isInitialized) AtomType::initialize();
    return AtomType::atomTypes[atomicNumber];
}

AtomType *AtomType::atomTypeFromAtomType(AtomTypes atomType) {
    return AtomType::atomTypeFromAtomicNumber(int(atomType));
}
string AtomType::name() const
{
    return m_name;
}

string AtomType::symbol() const
{
    return m_symbol;
}

double AtomType::massInverse() const
{
    return m_massInverse;
}

void AtomType::initialize()
{
    AtomType::atomTypes.push_back(new AtomType(0,1, "ERROR - ATOM NOT IN USE", "ER"));
    AtomType::atomTypes.push_back(new AtomType(1,1.00794, "Hydrogen", "H"));
    AtomType::atomTypes.push_back(new AtomType(2,4.002602, "Helium", "He"));
    AtomType::atomTypes.push_back(new AtomType(3,6.941, "Lithium", "Li"));
    AtomType::atomTypes.push_back(new AtomType(4,9.012182, "Beryllium", "Be"));
    AtomType::atomTypes.push_back(new AtomType(5,10.811, "Boron", "B"));
    AtomType::atomTypes.push_back(new AtomType(6,12.0107, "Carbon", "C"));
    AtomType::atomTypes.push_back(new AtomType(7,14.0067, "Nitrogen", "N"));
    AtomType::atomTypes.push_back(new AtomType(8,15.9994, "Oxygen", "O"));
    AtomType::atomTypes.push_back(new AtomType(9,18.9984032, "Fluorine", "F"));
    AtomType::atomTypes.push_back(new AtomType(10,20.1797, "Neon", "Ne"));
    AtomType::atomTypes.push_back(new AtomType(11,22.98976928, "Sodium", "Na"));
    AtomType::atomTypes.push_back(new AtomType(12,24.3050, "Magnesium", "Mg"));
    AtomType::atomTypes.push_back(new AtomType(13,26.9815386, "Aluminium", "Al"));
    AtomType::atomTypes.push_back(new AtomType(14,28.0855, "Silicon", "Si"));
    AtomType::atomTypes.push_back(new AtomType(15,30.973762, "Phosphorus", "P"));
    AtomType::atomTypes.push_back(new AtomType(16,32.065, "Sulfur", "S"));
    AtomType::atomTypes.push_back(new AtomType(17,35.453, "Chlorine", "Cl"));
    AtomType::atomTypes.push_back(new AtomType(18,39.948, "Argon", "Ar"));
}

int AtomType::atomicNumber() const
{
    return m_atomicNumber;
}

double AtomType::mass() const
{
    return m_mass;
}
