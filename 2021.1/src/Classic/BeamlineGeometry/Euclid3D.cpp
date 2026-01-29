// ------------------------------------------------------------------------
// $RCSfile: Euclid3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Euclid3D
//   Euclid3D represents an arbitrary 3-dimension translation and
//   rotation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Euclid3D.h"
#include <cmath>


// Class Euclid3D
// ------------------------------------------------------------------------
// The present implementation contains a Vector3D representing the
// displacement, and a Matrix3D representing the rotation.

Euclid3D::Euclid3D
(double x, double y, double z, double vx, double vy, double vz):
    V(x, y, z), R(vx, vy, vz), is_identity(false)
{}


Euclid3D::Euclid3D(const Vector3D &vec, const Rotation3D &rot):
    V(vec), R(rot), is_identity(false)
{}


bool Euclid3D::operator==(const Euclid3D &rhs) const {
    return (V == rhs.V)  && (R == rhs.R);
}


bool Euclid3D::operator!=(const Euclid3D &rhs) const {
    return (V != rhs.V)  || (R != rhs.R);
}


void Euclid3D::getAll(double &x, double &y, double &z,
                      double &vx, double &vy, double &vz) const {
    V.getComponents(x, y, z);
    R.getAxis(vx, vy, vz);
}


const Vector3D &Euclid3D::getVector() const {
    return V;
}


const Rotation3D &Euclid3D::getRotation() const {
    return R;
}


void Euclid3D::setX(double x) {
    V.setX(x);
}


void Euclid3D::setY(double y) {
    V.setY(y);
}


void Euclid3D::setZ(double z) {
    V.setZ(z);
}


void Euclid3D::setDisplacement(const Vector3D &v) {
    if(! v.isZero()) is_identity = false;
    V = v;
}


void Euclid3D::setRotation(const Rotation3D &rot) {
    if(! rot.isIdentity()) is_identity = false;
    R = rot;
}


Euclid3D Euclid3D::dot(const Euclid3D &rhs) const {
    Euclid3D result(*this);
    return result.dotBy(rhs);
}


const Euclid3D &Euclid3D::dotBy(const Euclid3D &rhs) {
    V += R * rhs.V;
    R *= rhs.R;
    if(! rhs.is_identity) is_identity = false;
    return *this;
}


Euclid3D Euclid3D::operator*(const Euclid3D &rhs) const {
    Euclid3D result(*this);
    return result.dotBy(rhs);
}


const Euclid3D &Euclid3D::operator*=(const Euclid3D &rhs) {
    return dotBy(rhs);
}


Euclid3D Euclid3D::inverse() const {
    Rotation3D rot = R.inverse();
    return Euclid3D(rot * (- V), rot);
}


bool Euclid3D::isPureTranslation() const {
    return R.isIdentity();
}


bool Euclid3D::isPureXRotation() const {
    return V.isZero() && R.isPureXRotation();
}


bool Euclid3D::isPureYRotation() const {
    return V.isZero() && R.isPureYRotation();
}


bool Euclid3D::isPureZRotation() const {
    return V.isZero() && R.isPureZRotation();
}


Euclid3D Euclid3D::identity() {
    return Euclid3D();
}


Euclid3D Euclid3D::translation(double x, double y, double z) {
    return Euclid3D(Vector3D(x, y, z), Rotation3D::Identity());
}


Euclid3D Euclid3D::XRotation(double angle) {
    Euclid3D result = Euclid3D(Vector3D(), Rotation3D::XRotation(angle));
    result.is_identity = angle == 0.0;
    return result;
}


Euclid3D Euclid3D::YRotation(double angle) {
    Euclid3D result = Euclid3D(Vector3D(), Rotation3D::YRotation(angle));
    result.is_identity = angle == 0.0;
    return result;
}


Euclid3D Euclid3D::ZRotation(double angle) {
    Euclid3D result = Euclid3D(Vector3D(), Rotation3D::ZRotation(angle));
    result.is_identity = angle == 0.0;
    return result;
}
