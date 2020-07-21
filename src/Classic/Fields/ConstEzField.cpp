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
#include "Fields/ConstEzField.h"


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
