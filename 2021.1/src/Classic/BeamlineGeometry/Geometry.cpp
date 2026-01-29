// ------------------------------------------------------------------------
// $RCSfile: BGeometryBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BGeometryBase
//    Pure virtual base class for all Beamline Geometries
//
// ------------------------------------------------------------------------
// Class category: BeamlineBGeometryBase
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class BGeometryBase.
// ------------------------------------------------------------------------

BGeometryBase::~BGeometryBase()
{}


void BGeometryBase::setElementLength(double)
{}


double BGeometryBase::getOrigin() const {
    return getArcLength() / 2.0;
}


double BGeometryBase::getEntrance() const {
    return - getOrigin();
}


double BGeometryBase::getExit() const {
    return getArcLength() - getOrigin();
}


Euclid3D BGeometryBase::getTotalTransform() const {
    return getTransform(getExit(), getEntrance());
}


Euclid3D BGeometryBase::getTransform(double s) const {
    return getTransform(0.0, s);
}


Euclid3D BGeometryBase::getEntranceFrame() const {
    return getTransform(0.0, getEntrance());
}


Euclid3D BGeometryBase::getExitFrame() const {
    return getTransform(0.0, getExit());
}


Euclid3D BGeometryBase::getEntrancePatch() const {
    return Euclid3D::identity();
}


Euclid3D BGeometryBase::getExitPatch() const {
    return Euclid3D::identity();
}
