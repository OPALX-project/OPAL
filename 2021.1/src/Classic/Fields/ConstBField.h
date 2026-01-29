#ifndef CLASSIC_ConstBField_HH
#define CLASSIC_ConstBField_HH

// ------------------------------------------------------------------------
// $RCSfile: ConstBField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstBField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/StaticMagneticField.h"


// Class ConstBField
// ------------------------------------------------------------------------
/// A homogenous magnetostatic field.
//  A static magnetic field independent of (x,y,z).

class ConstBField: public StaticMagneticField {

public:

    /// Default constructor.
    //  Constructs null field.
    ConstBField();

    virtual ~ConstBField();

    /// Get component.
    //  Return the x-component of the magnetic field in T.
    virtual double getBx() const;

    /// Get component.
    //  Return the y-component of the magnetic field in T.
    virtual double getBy() const;

    /// Get component.
    //  Return the z-component of the magnetic field in T.
    virtual double getBz() const;

    /// Set component.
    //  Assign the x-component of the magnetic field in T.
    virtual void setBx(double);

    /// Set component.
    //  Assign the y-component of the magnetic field in T.
    virtual void setBy(double);

    /// Set component.
    //  Assign the z-component of the magnetic field in T.
    virtual void setBz(double);
};

#endif // CLASSIC_ConstBField_HH
