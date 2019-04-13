#ifndef MYLAR_H
#define MYLAR_H

#include "Physics/Material.h"

namespace Physics {
    class Mylar: public Material {
    public:
        Mylar():
            Material(6.702,
                     12.88,
                     1.4,
                     39.95,
                     78.7,
                     {{3.350, 1683, 1900, 2.513e-02}})
        { }
    };
}
#endif