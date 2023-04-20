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

#include "PyOpal/PyCore/PyOpalObject.h"

namespace PyOpal {
namespace PyOpalObjectNS {

std::map<AttributeType, std::string> attributeName = std::map<AttributeType, std::string>({
    {DOUBLE, "float"},
    {STRING, "string"},
    {PREDEFINED_STRING, "predefined string"},
    {UPPER_CASE_STRING, "upper case string"},
    {BOOL, "bool"},
    {INT, "int"},
    {FLOAT_LIST, "list of floats"}
});
} // PyOpalObjectNS
} // PyOpal
