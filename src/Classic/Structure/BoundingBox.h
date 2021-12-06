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

    static std::shared_ptr<BoundingBox> getBoundingBox(const std::vector<Vector_t>& positions);

    void enlargeToContainPosition(const Vector_t& position);
    void enlargeToContainBoundingBox(const BoundingBox& boundingBox);

    boost::optional<Vector_t> getIntersectionPoint(const Vector_t& position,
                                                   const Vector_t& direction) const;

    bool isInside(const Vector_t& position) const;
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

#endif