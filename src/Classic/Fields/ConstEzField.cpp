// ------------------------------------------------------------------------
// $RCSfile: ConstEzField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstEzField
//   A static electric field of constant value in z-direction.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstEzField.h"


// Class ConstEzField
// ------------------------------------------------------------------------

ConstEzField::ConstEzField()
{}


ConstEzField::~ConstEzField()
{}


EVector ConstEzField::Efield(const Point3D &) const {
    return EVector(0.0, 0.0, Ez);
}


EVector ConstEzField::Efield(const Point3D &/*P*/, double) const {
    return EVector(0.0, 0.0, Ez);
}


double ConstEzField::getEz() const {
    return Ez;
}


void ConstEzField::setEz(double value) {
    Ez = value;
}


void ConstEzField::scale(double scalar) {
    Ez *= scalar;
}
