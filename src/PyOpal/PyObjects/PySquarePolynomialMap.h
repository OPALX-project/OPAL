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

#ifndef PYOPAL_PySquarePolynomialMap_H
#define PYOPAL_PySquarePolynomialMap_H

#include <Python.h>

namespace interpolation {

// note following are in interpolation namespace
class SquarePolynomialVector;
}

namespace PySquarePolynomialMap {

/** PyPolynomialMap is the python implementation of the C++ PolynomialMap class

 Provides a multivariate polynomial object
 */
typedef struct {
    PyObject_HEAD;
    interpolation::SquarePolynomialVector* map;
} PyPolynomialMap;

/** _alloc allocates memory for PyPolynomialMap
 *
 *  @param type - pointer to aPyPolynomialMapType object, as defined in
 *         PyPolynomialMap.cc
 *
 *  returns a PyPolynomialMap* (cast as a PyObject*); caller owns this memory
 */
static PyObject *_alloc(PyTypeObject *type, Py_ssize_t nitems);

/** _init initialises an allocated PyPolynomialMap object
 *
 *  @param self an initialised PyPolynomialMap* cast as a PyObject*; caller owns
 *         this memory
 *  @param args not used
 *  @param kwds not used
 *
 *  @returns 0 on success; -1 on failure
 */
static int _init(PyObject* self, PyObject *args, PyObject *kwds);

/** deallocate memory
 *
 *  @params self an initialised PyPolynomialMap*; memory will be freed by this 
 *          function
 */
static void _dealloc(PyPolynomialMap * self);

/** synonym for dealloc */
static void _free(PyPolynomialMap * self);

/** Initialise polynomial_map module
 *
 *  This is called by import polynomial_map; it initialises the PolynomialMap type
 *  allowing user to construct and call methods on PolynomialMap objects
 */
PyMODINIT_FUNC initsquare_polynomial_map(void);

/** Get the PolynomialMap coefficients
 *
 *  \param self a PyPolynomialMap
 *  \param args arguments to the function (not used)
 *  \param kwds keyword arguments to the function - see docstring
 *
 *  \return a list of lists of floats; each corresponding to a value in the
 *          polynomial map
 */
static PyObject* get_coefficients_as_matrix(PyObject *self, PyObject *args,
                                                                PyObject *kwds);

/** Calculate a polynomial vector
 *
 *  \param self a PyPolynomialMap
 *  \param args arguments to the function (not used)
 *  \param kwds keyword arguments to the function - see docstring
 *
 *  \return a PyList of floats with ValueDimension
 */
static PyObject* evaluate(PyObject *self, PyObject *args,
                                                                PyObject *kwds);

/** Get the point dimension
 *
 *  \param self a PyPolynomialMap
 *  \param args arguments to the function (not used)
 *  \param kwds keyword arguments to the function - see docstring
 *
 *  \return a PyLong
 */
// static PyObject* point_dimension(PyObject *self, PyObject *args,
//                                                              PyObject *kwds);

/** Get the value dimension
 *
 *  \param self a PyPolynomialMap
 *  \param args arguments to the function (not used)
 *  \param kwds keyword arguments to the function - see docstring
 *
 *  \return a PyLong
 */
// static PyObject* value_dimension(PyObject *self, PyObject *args,
//                                                              PyObject *kwds);

/** Get the string representation of the PolynomialMap
 *
 *  \param self a PyPolynomialMap
 *  \param args arguments to the function (not used)
 *  \param kwds keyword arguments to the function - see docstring
 *
 *  \return a PyString
 */
// static PyObject* str(PyObject *self, PyObject *args, PyObject *kwds);

}  // namespace PySquarePolynomialMap

#endif  // PYOPAL_PySquarePolynomialMap_H
