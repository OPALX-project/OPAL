#ifndef BORONCARBIDE_H
#define BORONCARBIDE_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class BoronCarbide: public Material {
    public:
        BoronCarbide():
            Material(26,
                     55.25,
                     2.48,
                     50.14,
                     84.7,
                     {{3.963, 6065.0, 1243.0, 7.782e-3}})
        { }
    };
}
#endif