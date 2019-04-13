#ifndef GRAPHITE_H
#define GRAPHITE_H

#include "Physics/Material.h"

namespace Physics {
    class Graphite: public Material {
    public:
        Graphite():
            Material(6,
                     12.0107,
                     2.210,
                     42.70,
                     78.0,
                     {{2.601, 1.701e3, 1.279e3, 1.638e-2}})
        { }
    };
}
#endif