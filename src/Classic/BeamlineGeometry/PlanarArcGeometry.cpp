// ------------------------------------------------------------------------
// $RCSfile: PlanarArcGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PlanarArcGeometry
//    Defines a Geometry which is a simple arc in the XZ plane.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Utilities/GeneralClassicException.h"

#include <cmath>
#include <algorithm>

// File scope function for calculating general transformations around arcs.
namespace {
    Euclid3D ArcTransform(double l, double h) {
        Euclid3D t;

        if(h) {
            double phi = h * l;
            t = Euclid3D::YRotation(- phi);
            t.setX((cos(phi) - 1.0) / h);
            t.setZ(sin(phi) / h);
        } else {
            t.setZ(l);
        }

        return t;
    }
}


PlanarArcGeometry::~PlanarArcGeometry()
{}


double PlanarArcGeometry::getArcLength() const {
    return len;
}


double PlanarArcGeometry::getElementLength() const {
    return len;
}


double PlanarArcGeometry::getCurvature() const {
    return h;
}


double PlanarArcGeometry::getBendAngle() const {
    return angle;
}


void PlanarArcGeometry::setBendAngle(double phi) {
    angle = phi;
    if(len != 0.0) h = angle / len;
}


void PlanarArcGeometry::setCurvature(double hh) {
    if(len != 0.0) {
        h = hh;
        angle = h * len;
    }
}


void PlanarArcGeometry::setElementLength(double l) {
    if (l < 0.0) {
        throw GeneralClassicException("PlanarArcGeometry::setElementLength",
                                      "The length of an element has to be positive");
    }
    len = l;
    if(len != 0.0) angle = h * len;
}


double PlanarArcGeometry::getOrigin() const {
    return len / 2.0;
}


double PlanarArcGeometry::getEntrance() const {
    return - len / 2.0;
}


double PlanarArcGeometry::getExit() const {
    return len / 2.0;
}


Euclid3D PlanarArcGeometry::getTransform(double fromS, double toS) const {
    return ArcTransform(toS - fromS, h);
}


Euclid3D PlanarArcGeometry::getTotalTransform() const {
    return ArcTransform(len, h);
}


Euclid3D PlanarArcGeometry::getTransform(double s) const {
    return ArcTransform(s, h);
}


Euclid3D PlanarArcGeometry::getEntranceFrame() const {
    return ArcTransform(- len / 2.0, h);
}


Euclid3D PlanarArcGeometry::getExitFrame() const {
    return ArcTransform(len / 2.0, h);
}
