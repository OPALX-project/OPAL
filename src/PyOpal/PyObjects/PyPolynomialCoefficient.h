//
// Python API for PolynomialCoefficient (part of the multidimensional polynomial fitting routines)
//
// Copyright (c) 2008-2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef PYOPAL_PyPolynomialCoefficient_H
#define PYOPAL_PyPolynomialCoefficient_H

#include <Python.h>


namespace interpolation {

// note following are in interpolation namespace
class PolynomialCoefficient;
}

namespace PyPolynomialCoefficient {

/** PyPolynomialMap is the python implementation of the C++ PolynomialMap class
 *
 *  Provides a multivariate polynomial object
 */
typedef struct {
    PyObject_HEAD;
    interpolation::PolynomialCoefficient* coeff;
} PyCoefficient;

}

#endif