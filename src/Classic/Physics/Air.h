#ifndef AIR_H
#define AIR_H

#include "Physics/Material.h"

namespace Physics {
    class Air: public Material {
    public:
        Air():
            Material(7,
                     14,
                     0.0012,
                     37.99,
                     85.7,
                     {{3.350, 1.683e3, 1.900e3, 2.513e-2}})
        { }
    };
}
#endif