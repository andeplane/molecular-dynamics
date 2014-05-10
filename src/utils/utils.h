#pragma once
#include <vector>
#include <iostream>

namespace CompPhys {
    namespace Utils {
        template <typename T>
        inline T &at(std::vector<T> &vec, int index) {
#ifdef DEBUG
            return vec.at(index);
#else
            return vec[index];
#endif
        }
    }
}
