#ifndef MSLANG_BOUNDINGBOX_H
#define MSLANG_BOUNDINGBOX_H

#include "Algorithms/Vektor.h"

#include <iostream>
#include <fstream>

namespace mslang {
    struct BoundingBox {
        Vector_t center_m;
        double width_m;
        double height_m;

        BoundingBox():
            center_m(0.0),
            width_m(0.0),
            height_m(0.0)
        { }

        BoundingBox(const BoundingBox &right):
            center_m(right.center_m),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        BoundingBox(const Vector_t &llc,
                    const Vector_t &urc):
            center_m(0.5 * (llc + urc)),
            width_m(urc[0] - llc[0]),
            height_m(urc[1] - llc[1])
        { }

        BoundingBox& operator=(const BoundingBox&) = default;
        bool doesIntersect(const BoundingBox &bb) const;
        bool isInside(const Vector_t &X) const;
        bool isInside(const BoundingBox &b) const;
        virtual void writeGnuplot(std::ostream &out) const;
        void print(std::ostream &out) const;
    };

    std::ostream & operator<< (std::ostream &out, const BoundingBox &bb);
}

#endif