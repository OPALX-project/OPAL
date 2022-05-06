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
#include "Utilities/MSLang/BoundingBox2D.h"

namespace mslang {
    bool BoundingBox2D::doesIntersect(const BoundingBox2D &bb) const {
        return ((center_m[0] - 0.5 * width_m <= bb.center_m[0] + 0.5 * bb.width_m) &&
                (center_m[0] + 0.5 * width_m >= bb.center_m[0] - 0.5 * bb.width_m) &&
                (center_m[1] - 0.5 * height_m <= bb.center_m[1] + 0.5 * bb.height_m) &&
                (center_m[1] + 0.5 * height_m >= bb.center_m[1] - 0.5 * bb.height_m));
    }

    bool BoundingBox2D::isInside(const Vector_t &X) const {
        if (2 * std::abs(X[0] - center_m[0]) <= width_m &&
            2 * std::abs(X[1] - center_m[1]) <= height_m)
            return true;

        return false;
    }

    bool BoundingBox2D::isInside(const BoundingBox2D &b) const {
        return (isInside(b.center_m + 0.5 * Vector_t( b.width_m,  b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t(-b.width_m,  b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t(-b.width_m, -b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t( b.width_m, -b.height_m, 0.0)));
    }

    void BoundingBox2D::writeGnuplot(std::ostream &out) const {
        std::vector<Vector_t> pts({Vector_t(center_m[0] + 0.5 * width_m, center_m[1] + 0.5 * height_m, 0),
                    Vector_t(center_m[0] - 0.5 * width_m, center_m[1] + 0.5 * height_m, 0),
                    Vector_t(center_m[0] - 0.5 * width_m, center_m[1] - 0.5 * height_m, 0),
                    Vector_t(center_m[0] + 0.5 * width_m, center_m[1] - 0.5 * height_m, 0)});
        unsigned int width = out.precision() + 8;
        for (unsigned int i = 0; i < 5; ++ i) {
            Vector_t & pt = pts[i % 4];

            out << std::setw(width) << pt[0]
                << std::setw(width) << pt[1]
                << std::endl;
        }
        out << std::endl;
    }

    void BoundingBox2D::print(std::ostream &out) const {
        out << std::setw(18) << center_m[0] - 0.5 * width_m
            << std::setw(18) << center_m[1] - 0.5 * height_m
            << std::setw(18) << center_m[0] + 0.5 * width_m
            << std::setw(18) << center_m[1] + 0.5 * height_m
            << std::endl;
    }

    std::ostream & operator<< (std::ostream &out, const BoundingBox2D &bb) {
        bb.print(out);
        return out;
    }
}