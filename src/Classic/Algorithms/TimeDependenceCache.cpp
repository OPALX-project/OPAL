//
// Class SinusoidalTimeDependence
//   A time dependence class that generates sine waves
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

#include "TimeDependenceCache.h"

void TimeDependenceCache::setTimeDependence(AbstractTimeDependence* timeDependence) {
    timeDependence_m = timeDependence;
    reset();
}

void TimeDependenceCache::reset() {
    head_m = 0;
    if (timeDependence_m) {
        // Start with the cache full to make expiry handling easier
        for (auto& [time_m, value_m, integral_m] : cache_m) {
            time_m     = 0.0;
            value_m    = timeDependence_m->getValue(0.0);
            integral_m = timeDependence_m->getIntegral(0.0);
        }
    }
}

std::array<TimeDependenceCache::Item, TimeDependenceCache::CacheSize>::iterator
TimeDependenceCache::find(double time) {
    return std::find_if(cache_m.begin(), cache_m.end(), [time](const Item& item) {
        return item.time_m == time;
    });
}

double TimeDependenceCache::getValue(const double time) {
    double result{};
    if (timeDependence_m) {
        if (const auto pos = find(time); pos == cache_m.end()) {
            result = placeInCache(time).value_m;
        } else {
            result = pos->value_m;
        }
    }
    return result;
}

double TimeDependenceCache::getIntegral(const double time) {
    double result{};
    if (timeDependence_m) {
        if (const auto pos = find(time); pos == cache_m.end()) {
            result = placeInCache(time).integral_m;
        } else {
            result = pos->integral_m;
        }
    }
    return result;
}

const TimeDependenceCache::Item& TimeDependenceCache::placeInCache(const double time) {
    auto& result      = cache_m.at(head_m);
    result.time_m     = time;
    result.value_m    = timeDependence_m->getValue(time);
    result.integral_m = timeDependence_m->getIntegral(time);
    ++head_m;
    if (head_m >= CacheSize) {
        head_m = 0;
    }
    return result;
}
