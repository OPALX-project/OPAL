#include "Utilities/MSLang/BoundingBox.h"

namespace mslang {
    bool BoundingBox::doesIntersect(const BoundingBox &bb) const {
        return ((center_m[0] - 0.5 * width_m <= bb.center_m[0] + 0.5 * bb.width_m) &&
                (center_m[0] + 0.5 * width_m >= bb.center_m[0] - 0.5 * bb.width_m) &&
                (center_m[1] - 0.5 * height_m <= bb.center_m[1] + 0.5 * bb.height_m) &&
                (center_m[1] + 0.5 * height_m >= bb.center_m[1] - 0.5 * bb.height_m));
    }

    bool BoundingBox::isInside(const Vector_t &X) const {
        if (2 * std::abs(X[0] - center_m[0]) <= width_m &&
            2 * std::abs(X[1] - center_m[1]) <= height_m)
            return true;

        return false;
    }

    bool BoundingBox::isInside(const BoundingBox &b) const {
        return (isInside(b.center_m + 0.5 * Vector_t( b.width_m,  b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t(-b.width_m,  b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t(-b.width_m, -b.height_m, 0.0)) &&
                isInside(b.center_m + 0.5 * Vector_t( b.width_m, -b.height_m, 0.0)));
    }

    void BoundingBox::writeGnuplot(std::ostream &out) const {
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

    void BoundingBox::print(std::ostream &out) const {
        out << std::setw(18) << center_m[0] - 0.5 * width_m
            << std::setw(18) << center_m[1] - 0.5 * height_m
            << std::setw(18) << center_m[0] + 0.5 * width_m
            << std::setw(18) << center_m[1] + 0.5 * height_m
            << std::endl;
    }

    std::ostream & operator<< (std::ostream &out, const BoundingBox &bb) {
        bb.print(out);
        return out;
    }
}