// ------------------------------------------------------------------------
// $RCSfile: BDipoleField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BDipoleField
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

#include "Fields/BDipoleField.h"


// Class BDipoleField
// ------------------------------------------------------------------------

BDipoleField::BDipoleField():
    Bx(0.0), By(0.0)
{}

BDipoleField::~BDipoleField()
{}


BVector BDipoleField::Bfield(const Point3D &) const {
    return BVector(Bx, By, 0.0);
}


BVector BDipoleField::Bfield(const Point3D &/*X*/, double) const {
    return BVector(Bx, By, 0.0);
}


double BDipoleField::getBx() const {
    return Bx;
}

double BDipoleField::getBy() const {
    return By;
}


void BDipoleField::setBx(double value) {
    Bx = value;
}

void BDipoleField::setBy(double value) {
    By = value;
}


BDipoleField &BDipoleField::addField(const BDipoleField &field) {
    Bx += field.Bx;
    By += field.By;
    return *this;
}


BDipoleField &BDipoleField::subtractField(const BDipoleField &field) {
    Bx -= field.Bx;
    By -= field.By;
    return *this;
}


void BDipoleField::scale(double scalar) {
    Bx *= scalar;
    By *= scalar;
}
