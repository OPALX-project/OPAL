#ifndef MOLYBDENUM_H
#define MOLYBDENUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Molybdenum: public Material {
    public:
        Molybdenum():
            Material(42,
                     95.94,
                     10.22,
                     9.8,
                     424.0,
                     {{7.248, 9.545e3, 4.802e2, 5.376e-3}})
        { }
    };
}
#endif