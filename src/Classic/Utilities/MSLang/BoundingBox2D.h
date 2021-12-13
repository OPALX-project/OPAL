//
// Class BoundingBox2D
//
// This class provides functionality to compute bounding boxes, to compute if a position
// is inside the box.
//
// Copyright (c) 2018 - 2021, Christof Metzger-Kraus
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef MSLANG_BOUNDINGBOX_H
#define MSLANG_BOUNDINGBOX_H

#include "Algorithms/Vektor.h"

#include <iostream>
#include <fstream>

namespace mslang {
    struct BoundingBox2D {
        Vector_t center_m;
        double width_m;
        double height_m;

        BoundingBox2D():
            center_m(0.0),
            width_m(0.0),
            height_m(0.0)
        { }

        BoundingBox2D(const BoundingBox2D &right):
            center_m(right.center_m),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        BoundingBox2D(const Vector_t &llc,
                    const Vector_t &urc):
            center_m(0.5 * (llc + urc)),
            width_m(urc[0] - llc[0]),
            height_m(urc[1] - llc[1])
        { }

        BoundingBox2D& operator=(const BoundingBox2D&) = default;
        bool doesIntersect(const BoundingBox2D &bb) const;
        bool isInside(const Vector_t &X) const;
        bool isInside(const BoundingBox2D &b) const;
        virtual void writeGnuplot(std::ostream &out) const;
        void print(std::ostream &out) const;
    };

    std::ostream & operator<< (std::ostream &out, const BoundingBox2D &bb);
}

#endif