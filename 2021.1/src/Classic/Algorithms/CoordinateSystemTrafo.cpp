#include "Algorithms/CoordinateSystemTrafo.h"
#include "Utility/Inform.h"

extern Inform *gmsg;

CoordinateSystemTrafo::CoordinateSystemTrafo():
    origin_m(0.0),
    orientation_m(1.0, 0.0, 0.0, 0.0),
    rotationMatrix_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{ }

CoordinateSystemTrafo::CoordinateSystemTrafo(const CoordinateSystemTrafo &right):
    origin_m(right.origin_m),
    orientation_m(right.orientation_m),
    rotationMatrix_m(right.rotationMatrix_m)
{ }

CoordinateSystemTrafo::CoordinateSystemTrafo(const Vector_t &origin,
                                             const Quaternion &orientation):
    origin_m(origin),
    orientation_m(orientation),
    rotationMatrix_m(orientation_m.getRotationMatrix())
{ }

void CoordinateSystemTrafo::invert() {
    origin_m = -orientation_m.rotate(origin_m);
    orientation_m = orientation_m.conjugate();
    rotationMatrix_m = transpose(rotationMatrix_m);
}

CoordinateSystemTrafo CoordinateSystemTrafo::operator*(const CoordinateSystemTrafo &right) const {
    CoordinateSystemTrafo result(*this);

    result *= right;
    return result;
}

void CoordinateSystemTrafo::operator*=(const CoordinateSystemTrafo &right) {
    origin_m = right.orientation_m.conjugate().rotate(origin_m) + right.origin_m;
    orientation_m *= right.orientation_m;
    orientation_m.normalize();
    rotationMatrix_m = orientation_m.getRotationMatrix();
}