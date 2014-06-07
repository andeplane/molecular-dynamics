#include <atomtype.h>
#include <iostream>

bool AtomType::isInitialized = false;
vector<shared_ptr<AtomType>> AtomType::atomTypes;

AtomType::AtomType(int atomicNumber, double mass, string name, string symbol) :
    m_atomicNumber(atomicNumber),
    m_mass(mass),
    m_massInverse(1.0/mass),
    m_name(name),
    m_symbol(symbol)
{

}

shared_ptr<AtomType> AtomType::atomTypeFromAtomicNumber(int atomicNumber) {
    if(!AtomType::isInitialized) AtomType::initialize();
    return AtomType::atomTypes[atomicNumber];
}

shared_ptr<AtomType> AtomType::atomTypeFromAtomType(AtomTypes atomType) {
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

void AtomType::initialize()
{

    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(0,1822.88841781, "ERROR - ATOM NOT IN USE", "ER")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(1,1837.36215184, "Hydrogen", "H")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(2,7296.29682689, "Helium", "He")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(3,12652.668508, "Lithium", "Li")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(4,16428.202187, "Beryllium", "Be")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(5,19707.2466849, "Boron", "B")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(6,21894.1659197, "Carbon", "C")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(7,25532.6512017, "Nitrogen", "N")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(8,29165.1209518, "Oxygen", "O")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(9,34631.9691501, "Fluorine", "F")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(10,36785.3414048, "Neon", "Ne")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(11,41907.7841485, "Sodium", "Na")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(12,44305.3029948, "Magnesium", "Mg")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(13,49184.3342085, "Aluminium", "Al")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(14,51196.7326583, "Silicon", "Si")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(15,56461.7120057, "Phosphorus", "P")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(16,58450.917117, "Sulfur", "S")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(17,64626.8630765, "Chlorine", "Cl")));
    AtomType::atomTypes.push_back(shared_ptr<AtomType>(new AtomType(18,72820.7465145, "Argon", "Ar")));
    AtomType::isInitialized = true;
}
