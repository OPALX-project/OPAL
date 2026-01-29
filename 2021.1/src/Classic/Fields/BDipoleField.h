#ifndef CLASSIC_BDipoleField_HH
#define CLASSIC_BDipoleField_HH

// ------------------------------------------------------------------------
// $RCSfile: BDipoleField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BDipoleField
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


// Class BDipoleField
// ------------------------------------------------------------------------
/// The field of a magnetic dipole.
//  A static magnetic dipole field in the (x,y)-plane.

class BDipoleField: public ConstBField {

public:

    /// Default constructor.
    //  Constructs a null field.
    BDipoleField();

    virtual ~BDipoleField();

    /// Get field.
    //  Return the time-independent part of the magnetic field in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P, double t) const;

    /// Get horizontal component.
    //  Return the horizontal component of the field in Teslas.
    virtual double getBx() const;

    /// Get vertical component.
    //  Return the vertical component of the field in Teslas.
    virtual double getBy() const;

    /// Set horizontal component.
    //  Assign the horizontal component of the field in Teslas.
    virtual void setBx(double Bx);

    /// Set vertical component.
    //  Assign the vertical component of the field in Teslas.
    virtual void setBy(double By);

    /// Add to field.
    //  Add [b]field[/b] to the old value; return new value.
    BDipoleField &addField(const BDipoleField &field);

    /// Subtract from field.
    //  Subtract [b]field[/b] from the old value; return new value.
    BDipoleField &subtractField(const BDipoleField &field);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    virtual void scale(double scalar);

private:

    // The field components.
    double Bx, By;
};

#endif // CLASSIC_BDipoleField_HH
