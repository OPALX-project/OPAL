#ifndef BERYLLIUM_H
#define BERYLLIUM_H

#include "Physics/Material.h"

namespace Physics {
    class Beryllium: public Material {
    public:
        Beryllium():
            Material(4,
                     9.012,
                     1.848,
                     65.19,
                     63.7,
                     {{2.590, 9.660e2, 1.538e2, 3.475e-2}})
        { }
    };
}
#endif