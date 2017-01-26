#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "core.hpp"

namespace mp {

    class vector : public container {
        double* data;
        unsigned int size;
    public:
        Vector() { }
        Vector(const Vector& v) {  }
    };

}

#endif
