//
// Class ValueRange
//
// This class provides functionality to compute the limits of a range of double values.
//
// Copyright (c) 2021, Christof Metzger-Kraus
//
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
#include "ValueRange.h"

#include <limits>
#include <algorithm>

ValueRange::ValueRange():
    minValue_m(std::numeric_limits<double>::max()),
    maxValue_m(std::numeric_limits<double>::lowest())
{}

void ValueRange::enlargeIfOutside(double value)
{
    minValue_m = std::min(minValue_m, value);
    maxValue_m = std::max(maxValue_m, value);
}