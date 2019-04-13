#ifndef GRAPHITER6710_H
#define GRAPHITER6710_H

#include "Physics/Material.h"

namespace Physics {
    class GraphiteR6710: public Material {
    public:
        GraphiteR6710():
            Material(6,
                     12.0107,
                     1.88,
                     42.70,
                     78.0,
                     {{2.601, 1.701e3, 1.279e3, 1.638e-2}})
        { }
    };
}
#endif