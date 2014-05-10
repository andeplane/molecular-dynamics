#pragma once
#include <iostream>
#include <sys/types.h>
#include <vector>
using std::cout;
using std::endl;

namespace CompPhys
{

template <class T>
class vector : public std::vector<T>
{
public:
    typedef typename std::vector<T>::const_reference const_reference;
    typedef typename std::vector<T>::reference reference;
    typedef typename std::vector<T>::size_type size_type;

    const_reference at(size_type __n) const
    {
#ifdef DEBUG
        if (__n >= this->size())
            throw std::out_of_range("CompPhys::vector out of range");
#endif
        return (*this)[__n];
    }

    reference at(size_type __n)
    {
#ifdef DEBUG
        if (__n >= this->size())
            throw std::out_of_range("CompPhys::vector out of range");
#endif
        return (*this)[__n];
    }
};

}
