#include "BoundingBox.h"

#include <limits>

BoundingBox::BoundingBox():
    lowerLeftCorner_m(std::numeric_limits<double>::max()),
    upperRightCorner_m(std::numeric_limits<double>::lowest())
{}

std::shared_ptr<BoundingBox> BoundingBox::getBoundingBox(const std::vector<Vector_t>& positions)
{
    std::shared_ptr<BoundingBox> boundingBox = std::make_shared<BoundingBox>();
    for (const Vector_t& position : positions) {
        boundingBox->enlargeToContainPosition(position);
    }

    return boundingBox;
}

void BoundingBox::enlargeToContainPosition(const Vector_t& position)
{
    for (unsigned int d = 0; d < 3; ++ d) {
        lowerLeftCorner_m[d] = std::min(lowerLeftCorner_m[d], position[d]);
        upperRightCorner_m[d] = std::max(upperRightCorner_m[d], position[d]);
    }
}

void BoundingBox::enlargeToContainBoundingBox(const BoundingBox& boundingBox)
{
    for (unsigned int d = 0; d < 3; ++ d) {
        lowerLeftCorner_m[d] = std::min(lowerLeftCorner_m[d], boundingBox.lowerLeftCorner_m[d]);
        upperRightCorner_m[d] = std::max(upperRightCorner_m[d], boundingBox.upperRightCorner_m[d]);
    }
}

boost::optional<Vector_t> BoundingBox::getIntersectionPoint(const Vector_t& position,
                                                            const Vector_t& direction) const
{
    boost::optional<Vector_t> result = boost::none;
    double minDistance = std::numeric_limits<double>::max();
    const Vector_t dimensions = upperRightCorner_m - lowerLeftCorner_m;
    Vector_t normal(1, 0, 0);
    for (unsigned int d = 0; d < 3; ++ d) {
        double sign = -1;
        Vector_t upperCorner = lowerLeftCorner_m + dot(normal, upperRightCorner_m) * normal;
        for (const Vector_t& p0 : {lowerLeftCorner_m, upperCorner}) {
            double tau = dot(p0 - position, sign * normal) / dot(direction, sign * normal);
            if (tau < 0.0) {
                continue;
            }
            Vector_t pointOnPlane = position + tau * direction;
            Vector_t relativeP = pointOnPlane - p0;
            bool isOnFace = true;
            for (unsigned int i = 1; i < 3; ++ i) {
                unsigned int idx = (d + i) % 3;
                if (relativeP[idx] < 0 ||
                    relativeP[idx] > dimensions[idx]) {
                    isOnFace = false;
                    break;
                }
            }
            if (isOnFace) {
                double distance = euclidean_norm(pointOnPlane - position);
                if (distance < minDistance) {
                    minDistance = distance;
                    result = pointOnPlane;
                }
            }
            sign *= -1;
        }

        normal = Vector_t(normal[2], normal[0], normal[1]);
    }

    return result;
}

bool BoundingBox::isInside(const Vector_t& position) const
{
    Vector_t relPosition = position - lowerLeftCorner_m;
    Vector_t dimensions = upperRightCorner_m - lowerLeftCorner_m;
    for (unsigned int d = 0; d < 3; ++ d) {
        if (relPosition[d] > dimensions[d]) return false;
    }
    return true;
}

void BoundingBox::print(std::ostream& output) const
{
    output << std::setprecision(8);
    output << "Bounding box\n"
           << "lower left corner: " << lowerLeftCorner_m << "\n"
           << "upper right corner: " << upperRightCorner_m << std::endl;
}