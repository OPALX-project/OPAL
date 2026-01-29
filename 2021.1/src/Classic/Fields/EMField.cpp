// ------------------------------------------------------------------------
// $RCSfile: EMField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Classes:
//   Point3D:      a point in 3-dimensional space.
//   EVector:      an electric field vector.
//   BVector:      a magnetic field vector.
//   EBVectors:    an electromagnetic field, containing E and B.
//   EMField:      a virtual base class for electromagnetic fields.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/EMField.h"


// Class Point3D
// ------------------------------------------------------------------------

Point3D::Point3D(double xx, double yy, double zz)
{ x = xx; y = yy; z = zz; }

double Point3D::getX() const
{ return x; }

double Point3D::getY() const
{ return y; }

double Point3D::getZ() const
{ return z; }


// Class EVector
// ------------------------------------------------------------------------

EVector::EVector(double Exx, double Eyy, double Ezz)
{ Ex = Exx; Ey = Eyy; Ez = Ezz; }


EVector EVector::operator*(double scale) const
{ return EVector(scale * Ex, scale * Ey, scale * Ez); }


double EVector::getEx() const
{ return Ex; }


double EVector::getEy() const
{ return Ey; }


double EVector::getEz() const
{ return Ez; }


// Class BVector
// ------------------------------------------------------------------------

BVector::BVector(double Bxx, double Byy, double Bzz)
{ Bx = Bxx; By = Byy; Bz = Bzz; }


BVector BVector::operator*(double scale) const
{ return BVector(scale * Bx, scale * By, scale * Bz); }


double BVector::getBx() const
{ return Bx; }


double BVector::getBy() const
{ return By; }


double BVector::getBz() const
{ return Bz; }


// Class EBVectors
// ------------------------------------------------------------------------

EBVectors::EBVectors(const EVector &eField, const BVector &bField):
    E(eField), B(bField)
{}


EVector EBVectors::getE() const
{ return E; }


double EBVectors::getEx() const
{ return E.getEx(); }


double EBVectors::getEy() const
{ return E.getEy(); }


double EBVectors::getEz() const
{ return E.getEz(); }


BVector EBVectors::getB() const
{ return B; }


double EBVectors::getBx() const
{ return B.getBx(); }


double EBVectors::getBy() const
{ return B.getBy(); }


double EBVectors::getBz() const
{ return B.getBy(); }


// Class EMField
// ------------------------------------------------------------------------

// Zero field representations.
const EVector
EMField::ZeroEfield(0.0, 0.0, 0.0);


const BVector
EMField::ZeroBfield(0.0, 0.0, 0.0);


const EBVectors
EMField::ZeroEBfield(EMField::ZeroEfield, EMField::ZeroBfield);


EMField::EMField()
{}


EMField::EMField(const EMField &)
{}


EMField::~EMField()
{}


const EMField &EMField::operator=(const EMField &) {
    return *this;
}


EVector EMField::Efield(const Point3D &) const {
    return ZeroEfield;
}


BVector EMField::Bfield(const Point3D &) const {
    return ZeroBfield;
}


EVector EMField::Efield(const Point3D &X, double) const {
    return Efield(X);
}


BVector EMField::Bfield(const Point3D &X, double) const {
    return Bfield(X);
}


EBVectors EMField::EBfield(const Point3D &X) const {
    return EBVectors(Efield(X), Bfield(X));
}


EBVectors EMField::EBfield(const Point3D &X, double) const {
    return EBfield(X);
}
