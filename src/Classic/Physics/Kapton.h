#ifndef KAPTON_H
#define KAPTON_H

#include "Physics/Material.h"

namespace Physics {
    class Kapton: public Material {
    public:
        Kapton():
            Material(6,
                     12,
                     1.4,
                     39.95,
                     79.6, // source: pstar from NIST
                     {{2.601, 1.701e3, 1.279e3, 1.638e-2}})
        { }
    };
}
#endif