#ifndef CLASSIC_BSingleMultipoleField_HH
#define CLASSIC_BSingleMultipoleField_HH

// ------------------------------------------------------------------------
// $RCSfile: BSingleMultipoleField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: BSingleMultipoleField
//   Representation for a single multipole field.
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
#include <cmath>
#include <complex>


// Template class BSingleMultipoleField
// ------------------------------------------------------------------------
/// Representation for a single magnetic multipole field.
//  The order and the skew flag are encoded in a template parameter.
//  {center}
//  order > 0: normal multipole, order < 0: skew multipole.
//  {/center}
//  Thus the compiler may optimize by unrolling loops on the order.
//  It should also omit some tests and eliminate unreachable code.

template <int order>
class BSingleMultipoleField: public EMField {

public:

    /// Default constructor.
    //  Constructs a null field.
    BSingleMultipoleField();

    BSingleMultipoleField(const BSingleMultipoleField &);
    virtual ~BSingleMultipoleField();
    BSingleMultipoleField &operator=(const BSingleMultipoleField &);

    /// Conversion operator.
    //  Return the field as a BMultipoleField.
    operator BMultipoleField() const;

    /// Field at a given point.
    //  Return the magnetic field at point [b]P[/b].
    virtual BVector Bfield(const Point3D &X) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    //  This default action returns the static part  BField(P).
    virtual BVector Bfield(const Point3D &P, double t) const;

    /// Return field coefficient.
    //  Return the single multipole coefficient in T/m**n.
    virtual double getComponent() const;

    /// Set field coefficient.
    //  Assign the single multipole coefficient in T/m**n.
    virtual void setComponent(double);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    void scale(double scalar);

    /// Return order.
    int size() const;

private:

    // The field component strength.
    double strength;
};


// Implementation
// ------------------------------------------------------------------------

template <int order> inline
BSingleMultipoleField<order>::BSingleMultipoleField():
    strength(0.0)
{}


template <int order> inline
BSingleMultipoleField<order>::BSingleMultipoleField
(const BSingleMultipoleField &field):
    strength(field.strength)
{}


template <int order>
BSingleMultipoleField<order>::~BSingleMultipoleField()
{}


template <int order> inline
BSingleMultipoleField<order> &BSingleMultipoleField<order>::
operator=(const BSingleMultipoleField &field) {
    strength = field.strength;
    return *this;
}


template <int order> inline
BSingleMultipoleField<order>::operator BMultipoleField() const {
    BMultipoleField field;
    if(order > 0) {
        field.setNormalComponent(order, strength);
    } else {
        field.setSkewComponent(order, strength);
    }
    return field;
}


template <int order> inline
BVector BSingleMultipoleField<order>::
Bfield(const Point3D &point) const {
    std::complex<double> z(point.getX(), point.getY());
    std::complex<double> B = strength;
    if(order > 0) {
        for(int i = 0; i < order; i++) {
            B *= z;
        }
        return BVector(std::real(B), -std::imag(B), 0.0);
    } else {
        for(int i = 0; i > order; i--) {
            B *= z;
        }
        return BVector(std::imag(B), std::real(B), 0.0);
    }
}


template <int order> inline
BVector BSingleMultipoleField<order>::
Bfield(const Point3D &point, double) const {
    return Bfield(point);
}


template <int order> inline
double BSingleMultipoleField<order>::getComponent() const {
    return strength;
}


template <int order> inline
void BSingleMultipoleField<order>::setComponent(double value) {
    strength = value;
}


template <int order> inline
void BSingleMultipoleField<order>::scale(double scalar) {
    strength *= scalar;
}


template <int order> inline
int BSingleMultipoleField<order>::size() const {
    return std::abs(order);
}

#endif // CLASSIC_BSingleMultipoleField_HH
