// ------------------------------------------------------------------------
// $RCSfile: ConstBzField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstBzField
//   A static magnetic field of constant value in z-direction.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstBzField.h"


// Class ConstBzField
// ------------------------------------------------------------------------

ConstBzField::ConstBzField()
{}


ConstBzField::~ConstBzField()
{}


BVector ConstBzField::Bfield(const Point3D &) const {
    return BVector(0.0, 0.0, Bz);
}


BVector ConstBzField::Bfield(const Point3D &/*X*/, double) const {
    return BVector(0.0, 0.0, Bz);
}


double ConstBzField::getBz() const {
    return Bz;
}

void ConstBzField::setBz(double value) {
    Bz = value;
}


void ConstBzField::scale(double scalar) {
    Bz *= scalar;
}
