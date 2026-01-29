// ------------------------------------------------------------------------
// $RCSfile: BMultipoleField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BMultipoleField
//   The static magnetic field of a multipole.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/BMultipoleField.h"
#include "Utilities/LogicalError.h"
#include <complex>


// Class BMultipoleField::Pair, inline functions.
// ------------------------------------------------------------------------
// These functions are used only within member functions of BMultipoleField
// and can be safely inlined.

inline BMultipoleField::Pair::Pair()
{ B = A = 0.0; }


inline BMultipoleField::Pair::Pair(double normal, double skewed)
{ B = normal; A = skewed; }


inline BMultipoleField::Pair::Pair(const Pair &rhs)
{ B = rhs.B; A = rhs.A; }


inline BMultipoleField::Pair::~Pair()
{}


inline void BMultipoleField::Pair::operator=(const Pair &rhs)
{ B = rhs.B; A = rhs.A; }


inline BMultipoleField::Pair
BMultipoleField::Pair::operator+(const Pair &rhs) const
{ return Pair(B + rhs.B, A + rhs.A); }


inline BMultipoleField::Pair
BMultipoleField::Pair::operator-(const Pair &rhs) const
{ return Pair(B - rhs.B, A - rhs.A); }


inline BMultipoleField::Pair
BMultipoleField::Pair::operator*(double scale) const
{ return Pair(B * scale, A * scale); }


inline void BMultipoleField::Pair::operator+=(const Pair &rhs)
{ B += rhs.B; A += rhs.A; }


inline void BMultipoleField::Pair::operator-=(const Pair &rhs)
{ B -= rhs.B; A -= rhs.A; }


inline void BMultipoleField::Pair::operator*=(double scale)
{ B *= scale; A *= scale; }


inline BMultipoleField::Pair BMultipoleField::Pair::operator-() const
{ return Pair(- B, - A); }


// Class BMultipoleField
// ------------------------------------------------------------------------

BMultipoleField::BMultipoleField():
    pairs(NULL),
    itsOrder(0)
{}


BMultipoleField::BMultipoleField(const BMultipoleField &rhs):
    pairs(NULL),
    itsOrder(rhs.itsOrder) {

    if (itsOrder > 0)
        pairs = new Pair[itsOrder];

    for(int i = 0; i < itsOrder; i++) {
        pairs[i] = rhs.pairs[i];
    }
}


BMultipoleField::~BMultipoleField() {
    delete [] pairs;
}


BMultipoleField &BMultipoleField::operator=(const BMultipoleField &rhs) {
    if(&rhs != this) {
        delete[] pairs;

        if (rhs.itsOrder > 0)
            pairs = new Pair[rhs.itsOrder];
        itsOrder = rhs.itsOrder;

        for(int i = 0; i < itsOrder; i++) {
            pairs[i] = rhs.pairs[i];
        }
    }

    return *this;
}


BVector BMultipoleField::Bfield(const Point3D &point) const {
    std::complex<double> X(point.getX(), point.getY());
    int i = itsOrder - 1;
    std::complex<double> B(pairs[i].B, pairs[i].A);

    while(i-- >= 0) {
        B = X * B / double(i + 1) + std::complex<double>(pairs[i].B, pairs[i].A);
    }

    return BVector(std::imag(B), std::real(B), 0.0);
}


BVector BMultipoleField::Bfield(const Point3D &point, double) const {
    return Bfield(point);
}


void BMultipoleField::setNormalComponent(int n, double b) {
    if(n <= 0) {
        throw LogicalError("BMultipoleField::setNormalComponent()",
                           "Field order should be > 0.");
    }

    if(n > itsOrder) reserve(n);
    pairs[n-1].B = b;
}


void BMultipoleField::setSkewComponent(int n, double a) {
    if(n <= 0) {
        throw LogicalError("BMultipoleField::setSkewComponent()",
                           "Field order should be > 0.");
    }

    if(n > itsOrder) reserve(n);
    pairs[n-1].A = a;
}


BMultipoleField &BMultipoleField::addField(const BMultipoleField &field) {
    if(field.itsOrder > itsOrder) reserve(field.itsOrder);
    for(int i = 0; i < field.itsOrder; ++i) pairs[i] += field.pairs[i];
    return *this;
}


BMultipoleField &BMultipoleField::subtractField(const BMultipoleField &field) {
    if(field.itsOrder > itsOrder) reserve(field.itsOrder);
    for(int i = 0; i < field.itsOrder; ++i) pairs[i] -= field.pairs[i];
    return *this;
}


void BMultipoleField::scale(double scalar) {
    for(int i = 0; i < itsOrder; ++i) {
        pairs[i] *= scalar;
    }
}


void BMultipoleField::reserve(int n) {
    if(n > itsOrder) {
        Pair *temp = new Pair[n];
        for(int i = 0; i < itsOrder; i++) temp[i] = pairs[i];
        delete [] pairs;
        itsOrder = n;
        pairs = temp;
    }
}
