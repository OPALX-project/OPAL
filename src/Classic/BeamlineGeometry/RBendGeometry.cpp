// ------------------------------------------------------------------------
// $RCSfile: RBendGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendGeometry
//    A Geometry which wraps a straight geometry to model the
//    geometry of an RBend element.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/RBendGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"
#include <cmath>


// Class RBendGeometry.
// ------------------------------------------------------------------------

RBendGeometry::RBendGeometry(double length, double angle):
    StraightGeometry(length), half_angle(angle / 2.0)
{}


RBendGeometry::RBendGeometry(const RBendGeometry &rhs):
    StraightGeometry(rhs), half_angle(rhs.half_angle)
{}


RBendGeometry::~RBendGeometry()
{}


double RBendGeometry::getArcLength() const {
    double length = StraightGeometry::getElementLength();
    return (half_angle == 0.0) ? length : length * half_angle / sin(half_angle);
}


double RBendGeometry::getBendAngle() const {
    return 2.0 * half_angle;
}


void RBendGeometry::setBendAngle(double angle) {
    half_angle = angle / 2.0;
}


Euclid3D RBendGeometry::getTotalTransform() const {
    Euclid3D patch = Euclid3D::YRotation(- half_angle);
    Euclid3D body  = StraightGeometry::getTotalTransform();
    return patch * body * patch;
}


Euclid3D RBendGeometry::getEntranceFrame() const {
    return
        StraightGeometry::getEntranceFrame() *
        Euclid3D::YRotation(half_angle);
}


Euclid3D RBendGeometry::getExitFrame() const {
    return
        StraightGeometry::getExitFrame() *
        Euclid3D::YRotation(- half_angle);
}


Euclid3D RBendGeometry::getEntrancePatch() const {
    double d = getElementLength() / 4.0 * tan(half_angle / 2.0);
    return
        Euclid3D::YRotation(- half_angle) * Euclid3D::translation(d, 0.0, 0.0);
}


Euclid3D RBendGeometry::getExitPatch() const {
    double d = getElementLength() / 4.0 * tan(half_angle / 2.0);
    return
        Euclid3D::translation(- d, 0.0, 0.0) * Euclid3D::YRotation(- half_angle);
}