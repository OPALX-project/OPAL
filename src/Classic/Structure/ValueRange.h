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

class ValueRange {
public:
    ValueRange();

    void enlargeIfOutside(double value);

    bool isInside(double value) const
    {
        return minValue_m < value && value < maxValue_m;
    }

    bool isOutside(double value) const
    {
        return !isInside(value);
    }

private:
    double minValue_m;
    double maxValue_m;
};

#endif