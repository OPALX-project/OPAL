#ifndef CLASSIC_ConstBzField_HH
#define CLASSIC_ConstBzField_HH

// ------------------------------------------------------------------------
// $RCSfile: ConstBzField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstBzField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstBField.h"


// Class ConstBzField
// ------------------------------------------------------------------------
/// A homogeneous magnetostatic field in z-direction.

class ConstBzField: public ConstBField {

public:

    /// Default constructor.
    //  Constructs a null field.
    ConstBzField();

    virtual ~ConstBzField();

    /// Get field.
    //  Return the time-independent part of the magnetic field in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P, double t) const;

    /// Get component.
    //  Return the z-component of the magnetic field in T.
    virtual double getBz() const;

    /// Set component.
    //  Assign the z-component of the magnetic field in T.
    virtual void setBz(double Bz);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    virtual void scale(double scalar);

private:

    // The field component.
    double Bz;
};

#endif // CLASSIC_ConstBzField_HH
