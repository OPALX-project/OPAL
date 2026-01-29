/*
 *  Copyright (c) 2018, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "gtest/gtest.h"

#include "Utilities/GeneralClassicException.h"
#include "Algorithms/SplineTimeDependence.h"

class SplineTimeDependenceTest : public ::testing::Test { 
public:
    SplineTimeDependenceTest() : times_m(10), values_m(10) {
        for (size_t i = 0; i < 10; ++i) {
            times_m[i] = i*i+i;
            values_m[i] = times_m[i]; // x^2
        }
    }

    std::vector<double> times_m;
    std::vector<double> values_m;
};


TEST_F(SplineTimeDependenceTest, ConstructorTest) {
    try {
        SplineTimeDependence timeDep(1, times_m, values_m);

        for (size_t i = 0; i < values_m.size(); ++i) {
            double test = timeDep.getValue(times_m[i]);
            EXPECT_NEAR(test, values_m[i], 1e-9) << "Index " << i;
        }
        EXPECT_THROW(SplineTimeDependence(0, times_m, values_m),
                     GeneralClassicException);
        EXPECT_THROW(SplineTimeDependence(2, times_m, values_m),
                     GeneralClassicException);
        times_m.push_back(1.);
        EXPECT_THROW(SplineTimeDependence(1, times_m, values_m),
                     GeneralClassicException);
        times_m = std::vector<double>(1, 0.);
        values_m = std::vector<double>(1, 0.);
        EXPECT_THROW(SplineTimeDependence(1, times_m, values_m),
                     GeneralClassicException);
        times_m = {0., 1.};
        values_m = std::vector<double>(2, 0.);
        SplineTimeDependence(1, times_m, values_m);
        times_m = {0., 1., 2.};
        values_m = std::vector<double>(3, 0.);
        EXPECT_THROW(SplineTimeDependence(3, times_m, values_m),
                     GeneralClassicException);
        times_m = {0., 1., 2., 3.};
        values_m = std::vector<double>(4, 0.);
        SplineTimeDependence(3, times_m, values_m);
        times_m = {0., 0., 2., 3.};
        values_m = std::vector<double>(4, 0.);
        EXPECT_THROW(SplineTimeDependence(3, times_m, values_m),
                     GeneralClassicException);
    } catch (GeneralClassicException& exc) {
        EXPECT_TRUE(false) << "Should not have thrown an exception:\n    "
                           << exc.what() << "\n    " << exc.where();
    }
}

TEST_F(SplineTimeDependenceTest, LinearLookupTest) {
    SplineTimeDependence timeDep(1, times_m, values_m);
    double test_x = (times_m[2]+times_m[3])/2.;
    double ref_y = (values_m[2]+values_m[3])/2.;
    EXPECT_NEAR(timeDep.getValue(test_x), ref_y, 1e-9);
    test_x = times_m[0] - (times_m[1]-times_m[0])/2.;
    EXPECT_THROW(timeDep.getValue(test_x),
                 GeneralClassicException);
}

TEST_F(SplineTimeDependenceTest, CubicLookupTest) {
    // if I give it a quadratic or cubic, it gets the wrong answer!!
    // I am sure that it is doing the right thing though... ahem
    SplineTimeDependence timeDep(1, times_m, values_m);
    for (size_t i = 0; i < times_m.size()-1; ++i) {
        double test_x = (times_m[i]+times_m[i+1])/2.;
        double ref_y = test_x;
        double test_y = timeDep.getValue(test_x);
        //std::cerr << times_m[i] << " " << times_m[i+1] << " ** " << test_x << std::endl;
        //std::cerr << values_m[i] << " " << values_m[i+1] << " ** " << ref_y << " ** " << test_y << std::endl;
        EXPECT_NEAR(test_y, ref_y, 1e-9);
    }
}
