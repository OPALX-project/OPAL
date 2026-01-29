// ------------------------------------------------------------------------
// $RCSfile: NullGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NullGeometry
//    A geometry representing an identity transform.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/NullGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class NullGeometry
// ------------------------------------------------------------------------

NullGeometry::~NullGeometry()
{}


double NullGeometry::getArcLength() const {
    return 0.0;
}


double NullGeometry::getElementLength() const {
    return 0.0;
}


double NullGeometry::getOrigin() const {
    return 0.0;
}


Euclid3D NullGeometry::getTransform(double, double) const {
    return Euclid3D::identity();
}


Euclid3D NullGeometry::getTransform(double) const {
    return Euclid3D::identity();
}

