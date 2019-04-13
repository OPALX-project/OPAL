#ifndef ALUMINAAL2O3_H
#define ALUMINAAL2O3_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class AluminaAL2O3: public Material {
    public:
        AluminaAL2O3():
            Material(50,
                     101.96,
                     3.97,
                     27.94,
                     145.2,
                     {{7.227, 1.121e4, 3.864e2, 4.474e-3}})
        { }
    };
}
#endif