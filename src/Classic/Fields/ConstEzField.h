//
// Class ConstEzField
//   A homogeneous electrostatic field in z-direction.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_ConstEzField_HH
#define CLASSIC_ConstEzField_HH

#include "Fields/StaticElectricField.h"

class ConstEzField: public StaticElectricField {

public:

    /// Default constructor.
    //  Constructs a null field.
    ConstEzField();

    virtual ~ConstEzField();

    /// Get field.
    //  Return the time-independent part of the electric field in point [b]P[/b].
    virtual EVector Efield(const Point3D &P) const;

    /// Get field.
    //  Return the electric field at time [b]t[/b] in point [b]P[/b].
    virtual EVector Efield(const Point3D &P, double t) const;

    /// Get component.
    //  Return the x-component of the electric field in A/m.
    virtual double getEz() const;

    /// Set component.
    //  Assign the z-component of the electric field in A/m.
    virtual void setEz(double);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    virtual void scale(double scalar);

private:

    // The field components.
    double Ez;
};

#endif // CLASSIC_ConstEzField_HH
