#ifndef TITANIUM_H
#define TITANIUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Titanium: public Material {
    public:
        Titanium():
            Material(22,
                     47.8,
                     4.54,
                     16.16,
                     233.0,
                     {{5.489, 5.260e3, 6.511e2, 8.930e-3}})
        { }
    };
}
#endif