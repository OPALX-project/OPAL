// ------------------------------------------------------------------------
// $RCSfile: EDipoleField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EDipoleField
//   An electrostatic dipole field in the (x,y)-plane.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/EDipoleField.h"


// Class EDipoleField
// ------------------------------------------------------------------------

EDipoleField::EDipoleField()
{}


EDipoleField::~EDipoleField()
{}


EVector EDipoleField::Efield(const Point3D &) const {
    return EVector(Ex, Ey, 0.0);
}


EVector EDipoleField::Efield(const Point3D &/*P*/, double) const {
    return EVector(Ex, Ey, 0.0);
}


double EDipoleField::getEx() const {
    return Ex;
}

double EDipoleField::getEy() const {
    return Ey;
}


void EDipoleField::setEx(double value) {
    Ex = value;
}


void EDipoleField::setEy(double value) {
    Ey = value;
}


void EDipoleField::scale(double scalar) {
    Ex *= scalar;
    Ey *= scalar;
}
