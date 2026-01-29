// ------------------------------------------------------------------------
// $RCSfile: Vector3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Vector3D
//   Vector3D represents an 3-dimension vector.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Vector3D.h"


// Class Vector3D
// ------------------------------------------------------------------------

Vector3D::Vector3D(double x, double y, double z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
}


Vector3D &Vector3D::operator+=(const Vector3D &rhs) {
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    return *this;
}


Vector3D &Vector3D::operator-=(const Vector3D &rhs) {
    v[0] -= rhs.v[0];
    v[1] -= rhs.v[1];
    v[2] -= rhs.v[2];
    return *this;
}


Vector3D &Vector3D::operator*=(double rhs) {
    v[0] *= rhs;
    v[1] *= rhs;
    v[2] *= rhs;
    return *this;
}


bool Vector3D::operator==(const Vector3D &rhs) const {
    return (v[0] == rhs.v[0]) && (v[1] == rhs.v[1]) && (v[0] == rhs.v[1]);
}


bool Vector3D::operator!=(const Vector3D &rhs) const {
    return (v[0] != rhs.v[0]) || (v[1] != rhs.v[1]) || (v[0] != rhs.v[1]);
}


Vector3D Vector3D::operator-() const {
    return Vector3D(- v[0], - v[1], - v[2]);
}


Vector3D operator+(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs(0) + rhs(0), lhs(1) + rhs(1), lhs(2) + rhs(2));
}


Vector3D operator-(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs(0) - rhs(0), lhs(1) - rhs(1), lhs(2) - rhs(2));
}


Vector3D operator*(const Vector3D &lhs, double rhs) {
    return Vector3D(lhs(0) * rhs, lhs(1) * rhs, lhs(2) * rhs);
}


Vector3D operator*(double lhs, const Vector3D &rhs) {
    return Vector3D(lhs * rhs(0), lhs * rhs(1), lhs * rhs(2));
}


void Vector3D::clear() {
    v[0] = v[1] = v[2] = 0.0;
}


void Vector3D::getComponents(double &x, double &y, double &z)
const {
    x = v[0];
    y = v[1];
    z = v[2];
}


bool Vector3D::isZero() const {
    return v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0;
}


Vector3D cross(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs(1) * rhs(2) - lhs(2) * rhs(1),
                    lhs(2) * rhs(0) - lhs(0) * rhs(2),
                    lhs(0) * rhs(1) - lhs(1) * rhs(0));
}


double dot(const Vector3D &lhs, const Vector3D &rhs) {
    return (lhs(0) * rhs(0) + lhs(1) * rhs(1) + lhs(2) * rhs(2));
}
