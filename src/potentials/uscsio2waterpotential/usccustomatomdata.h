#pragma once
#include <memory>

class USCCustomAtomData : public CustomAtomData
{
private:
    Atom *m_atom;
    double m_nHydrogen;
    double m_nSilicon;
public:
    USCCustomAtomData();
    virtual void reset();
};
