#pragma once

#include <vector>
#include <string>
using std::string; using std::vector;

enum class AtomTypes {NoAtom=0, Hydrogen = 1, Helium = 2, Lithium = 3, Beryllium = 4, Boron = 5, Carbon = 6, Nitrogen = 7, Oxygen = 8, Fluorine = 9, Neon = 10, Sodium = 11, Magnesium = 12, Aluminium = 13, Silicon = 14, Phosphorus = 15, Sulfur = 16, Chlorine = 17, Argon = 18};
inline int operator + (AtomTypes val) {
    return static_cast<int>(val);
}

class AtomType
{
private:
    static vector<AtomType*> atomTypes;
    static bool isInitialized;
    static void initialize();

public:
    static AtomType *atomTypeFromAtomicNumber(int atomicNumber);
    static AtomType *atomTypeFromAtomType(AtomTypes atomType);

    AtomType(int atomicNumber, double mass, string name, string symbol);

    inline int atomicNumber() const {
        return m_atomicNumber;
    }
    inline double mass() const {
        return m_mass;
    }
    inline double massInverse() const {
        return m_massInverse;
    }

    string name() const;
    string symbol() const;
private:
    int m_atomicNumber;
    double m_mass;
    double m_massInverse;
    string m_name;
    string m_symbol;
};
