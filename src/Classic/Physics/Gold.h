#ifndef GOLD_H
#define GOLD_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Gold: public Material {
    public:
        Gold():
            Material(79,
                     197,
                     19.3,
                     6.46,
                     790.0,
                     {{5.458, 7.852e3, 9.758e2, 2.077e-2}})
        { }
    };
}
#endif