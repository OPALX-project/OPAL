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

#ifndef PYOPAL_PYCORE_EXCEPTIONTRANSLATION_H
#define PYOPAL_PYCORE_EXCEPTIONTRANSLATION_H 1

#include <Python.h>
#include <structmember.h>

#include <string>
#include <exception>

#include <boost/python.hpp>

namespace PyOpal {
namespace py = boost::python;

/** Exception translation uses boost::python hooks to wrap C++ exceptions */
namespace ExceptionTranslation {

/** Register exception translations with boost */
void registerExceptions();

/** Translates the std::exceptions into a RuntimeError */
template <class T>
void translateException(T const& exception) {
    PyErr_SetString(PyExc_RuntimeError, exception.what());
}

/** Translates the OpalExceptions into a RuntimeError */
template <class T>
void translateOpalException(T const& exception) {
    std::string msg = exception.what()+" in C++ method "+exception.where();
    PyErr_SetString(PyExc_RuntimeError, msg.c_str());
} 


}
}

#endif
