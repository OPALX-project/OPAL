// ------------------------------------------------------------------------
// $RCSfile: StraightGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StraightGeometry
//    A concrete geometry representing a straight line segment which runs
//    coaxial with the local z-axis.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/StraightGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Utilities/GeneralClassicException.h"

#include <algorithm>

// Class StraightGeometry.
// ------------------------------------------------------------------------

StraightGeometry::~StraightGeometry()
{}


double StraightGeometry::getArcLength() const {
    return len;
}


double StraightGeometry::getElementLength() const {
    return len;
}


void StraightGeometry::setElementLength(double l) {
    if (l < 0.0) {
        throw GeneralClassicException("StraightGeometry::setElementLength",
                                      "The length of an element has to be positive");
    }
    len = std::max(0.0, l);
}


double StraightGeometry::getOrigin() const {
    return len / 2.0;
}


double StraightGeometry::getEntrance() const {
    return -len / 2.0;
}


double StraightGeometry::getExit() const {
    return len / 2.0;
}


Euclid3D StraightGeometry::getTransform(double fromS, double toS) const {
    return Euclid3D::translation(0, 0, fromS - toS);
}


Euclid3D StraightGeometry::getTotalTransform() const {
    return Euclid3D::translation(0, 0, len);
}


Euclid3D StraightGeometry::getTransform(double s) const {
    return Euclid3D::translation(0, 0, s);
}


Euclid3D StraightGeometry::getEntranceFrame() const {
    return Euclid3D::translation(0, 0, -len / 2.0);
}


Euclid3D StraightGeometry::getExitFrame() const {
    return Euclid3D::translation(0, 0, len / 2.0);
}