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

#include "Classic/Utilities/ClassicException.h"
#include "Utilities/OpalException.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"

namespace PyOpal {
namespace ExceptionTranslation {
void registerExceptions() {
    // handle std::exception (e.g. for IO errors)
    py::register_exception_translator<std::exception>(translateException<std::exception>);
    // handle Opal exceptions (they all inherit from ClassicException)
    py::register_exception_translator<ClassicException>(translateOpalException<ClassicException>);
}
}
}
