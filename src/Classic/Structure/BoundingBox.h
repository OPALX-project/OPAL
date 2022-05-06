//
// Class BoundingBox: defines a 3D box and provides functionality to determine whether a point
//                    is in- or outside it. Additionally provides functionality to compute the
//                    intersection point with a line.
//
// This class provides functionality to compute bounding boxes, to compute if a position
// is inside the box and to compute the intersection point between a ray and the bounding box.
//
// Copyright (c) 201x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include "Algorithms/Vektor.h"

#include "boost/optional.hpp"

#include <vector>
#include <utility>

class BoundingBox
{
public:
    BoundingBox();

    static BoundingBox getBoundingBox(const std::vector<Vector_t>& positions);

    void enlargeToContainPosition(const Vector_t& position);
    void enlargeToContainBoundingBox(const BoundingBox& boundingBox);

    boost::optional<Vector_t> getIntersectionPoint(const Vector_t& position,
                                                   const Vector_t& direction) const;

    bool isInside(const Vector_t& position) const;
    bool isOutside(const Vector_t& position) const;
    void print(std::ostream& output) const;

    std::pair<Vector_t, Vector_t> getCorners() const;
private:
    Vector_t lowerLeftCorner_m;
    Vector_t upperRightCorner_m;
};

inline
std::pair<Vector_t, Vector_t> BoundingBox::getCorners() const {
    return std::make_pair(lowerLeftCorner_m, upperRightCorner_m);
}

inline
bool BoundingBox::isOutside(const Vector_t& position) const {
    return !isInside(position);
}
#endif
