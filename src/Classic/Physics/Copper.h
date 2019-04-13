#ifndef COPPER_H
#define COPPER_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Copper: public Material {
    public:
        Copper():
            Material(29,
                     63.54,
                     8.96,
                     12.86,
                     322.0,
                     {{4.194, 4.649e3, 8.113e1, 2.242e-2}})
        { }
    };
}
#endif