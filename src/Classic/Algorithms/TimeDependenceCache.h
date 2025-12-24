//
// Class TimeDependenceCache
//   A class that accesses a time dependence through a short cache
//   The cache should be just long enough so that the Runge-Kutta calls
//   for a single step fit inside.
//   The cache consists of an array of values used in a circular buffer
//   manner, so the last few calls for a time are stored.
//   The cache is short to keep search times very low.
//
// Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
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

#ifndef OPAL_TIMEDEPENDENCECACHE_H
#define OPAL_TIMEDEPENDENCECACHE_H

#include "AbstractTimeDependence.h"
#include <array>

class TimeDependenceCache {
public:
    static constexpr size_t CacheSize = 6;

    struct Item {
        double time_m;
        double value_m;
        double integral_m;
    };

    void setTimeDependence(AbstractTimeDependence* timeDependence);
    double getValue(double time);
    double getIntegral(double time);
    void reset();
    size_t getHead() const { return head_m; }
    std::array<Item, CacheSize>::iterator find(double time);
    std::array<Item, CacheSize>::iterator begin() { return cache_m.begin(); }
    std::array<Item, CacheSize>::iterator end() { return cache_m.end(); }

private:
    AbstractTimeDependence* timeDependence_m{};

    std::array<Item, CacheSize> cache_m{};
    size_t head_m{};

    const Item& placeInCache(double time);
};

#endif  // OPAL_TIMEDEPENDENCECACHE_H
