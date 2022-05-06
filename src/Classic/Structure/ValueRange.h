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
#ifndef VALUERANGE_H
#define VALUERANGE_H

#include "Utility/Inform.h"

#include <algorithm>
#include <limits>

template<class T>
class ValueRange {
public:
    ValueRange():
        minValue_m(std::numeric_limits<T>::max()),
        maxValue_m(std::numeric_limits<T>::lowest())
    {}

    void enlargeIfOutside(T value)
    {
        minValue_m = std::min(minValue_m, value);
        maxValue_m = std::max(maxValue_m, value);
    }
    bool isInside(T value) const
    {
        return minValue_m < value && value < maxValue_m;
    }

    bool isOutside(T value) const
    {
        return !isInside(value);
    }

    void print(Inform& out) const
    {
        out << "Value range between " << minValue_m << " and " << maxValue_m;
    }
private:
    T minValue_m;
    T maxValue_m;
};

template<class T>
Inform& operator<<(Inform& out, ValueRange<T> const& range)
{
    range.print(out);
    return out;
}

#endif
