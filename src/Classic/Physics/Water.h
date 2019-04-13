#ifndef WATER_H
#define WATER_H

#include "Physics/Material.h"

namespace Physics {
    class Water: public Material {
    public:
        Water():
            Material(10,
                     18,
                     1,
                     36.08,
                     75.0, // source: pstar from NIST
                     {{2.199, 2.393e3, 2.699e3, 1.568e-2}})
        { }
    };
}
#endif