#ifndef ALUMINUM_H
#define ALUMINUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Aluminum: public Material {
    public:
        Aluminum():
            Material(13,
                     26.98,
                     2.7,
                     24.01,
                     166.0,
                     {{4.739, 2.766e3, 1.645e2, 2.023e-2}})
        { }
    };
}
#endif