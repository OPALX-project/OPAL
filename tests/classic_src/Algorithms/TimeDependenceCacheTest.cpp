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

#include "gtest/gtest.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Algorithms/TimeDependenceCache.h"
#include "opal_test_utilities/SilenceTest.h"


TEST(TimeDependenceCacheTest, Operation) {
    OpalTestUtilities::SilenceTest silencer;
    PolynomialTimeDependence polynomialTimeDependence{{1.0, 2.0, 3.0}};
    TimeDependenceCache cache;
    cache.setTimeDependence(&polynomialTimeDependence);
    // Time = 0.1 should not be in the cache
    EXPECT_EQ(cache.find(0.1), cache.end());
    // Check the result of getting the value
    EXPECT_DOUBLE_EQ(cache.getValue(0.1), 1.23);
    // Now 0.1 should be in the cache
    EXPECT_NE(cache.find(0.1), cache.end());
    EXPECT_DOUBLE_EQ(cache.getValue(0.1), 1.23);
    // Now get time = 0.2 which shouldn't be in the cache
    EXPECT_EQ(cache.find(0.2), cache.end());
    EXPECT_DOUBLE_EQ(cache.getValue(0.2), 1.52);
    // Now they should both be in the cache
    EXPECT_NE(cache.find(0.1), cache.end());
    EXPECT_NE(cache.find(0.2), cache.end());
    EXPECT_DOUBLE_EQ(cache.getValue(0.1), 1.23);
    EXPECT_DOUBLE_EQ(cache.getValue(0.2), 1.52);
    // Put 5 more values in the cache
    cache.getValue(0.3);
    cache.getValue(0.4);
    cache.getValue(0.5);
    cache.getValue(0.6);
    cache.getValue(0.7);
    // Which should mean the first one is no longer in the cache, but the second one is still there
    EXPECT_EQ(cache.find(0.1), cache.end());
    EXPECT_NE(cache.find(0.2), cache.end());
    // Resetting the cache should clear them both
    cache.reset();
    EXPECT_EQ(cache.find(0.1), cache.end());
    EXPECT_EQ(cache.find(0.2), cache.end());
}
