/*
 *  Copyright (c) 2025, Jon Thompson
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
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/SinusoidalTimeDependence.h"

#include "opal_test_utilities/SilenceTest.h"

TEST(SinusoidalTimeDependenceTest, SinusoidalTimeDependenceTest) {
    OpalTestUtilities::SilenceTest silencer;

    // Check empty coefficients always returns 0.
    SinusoidalTimeDependence time_dependence_1({}, {}, {}, {});
    EXPECT_DOUBLE_EQ(time_dependence_1.getValue(0.1), 0.0);

    // Check a sine wave value
    SinusoidalTimeDependence time_dependence_2({8.0}, {}, {}, {});
    EXPECT_DOUBLE_EQ(time_dependence_2.getValue(0.1), -0.47552825814757682);

    // Check the amplitude
    SinusoidalTimeDependence time_dependence_3({8.0}, {}, {2.0}, {});
    EXPECT_DOUBLE_EQ(time_dependence_3.getValue(0.1), -0.95105651629515364);

    // Check the phase offset
    SinusoidalTimeDependence time_dependence_4({8.0}, {0.1}, {2.0} , {});
    EXPECT_DOUBLE_EQ(time_dependence_4.getValue(0.1), -0.91545497277810161);

    // Check the DC offset
    SinusoidalTimeDependence time_dependence_5({8.0}, {0.1}, {2.0}, {-1.0});
    EXPECT_DOUBLE_EQ(time_dependence_5.getValue(0.1), -1.91545497277810161);

    // Check clone produces same result
    SinusoidalTimeDependence* time_dependence_clone = time_dependence_5.clone();
    EXPECT_DOUBLE_EQ(time_dependence_clone->getValue(0.1),-1.91545497277810161);
    delete time_dependence_clone;
}

TEST(SinusoidalTimeDependenceTest, TDMapTest) {
    OpalTestUtilities::SilenceTest silencer;

    // throw on empty value
    EXPECT_THROW(AbstractTimeDependence::getTimeDependence("name"),
                 GeneralClassicException);

    // set/get time dependence
    SinusoidalTimeDependence time_dep({}, {}, {}, {});
    std::shared_ptr<SinusoidalTimeDependence> td1(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    std::shared_ptr<SinusoidalTimeDependence> td2(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td2", td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td2"), td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    // set time dependence overwriting existing time dependence
    // should overwrite, without memory leak
    std::shared_ptr<SinusoidalTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td3);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td3);
}

TEST(SinusoidalTimeDependenceTest, TDMapNameLookupTest) {
    OpalTestUtilities::SilenceTest silencer;

    EXPECT_THROW(AbstractTimeDependence::getName(nullptr),
                 GeneralClassicException);
    SinusoidalTimeDependence time_dep({}, {}, {}, {});
    std::shared_ptr<SinusoidalTimeDependence> td1(time_dep.clone());
    std::shared_ptr<SinusoidalTimeDependence> td2(time_dep.clone());
    std::shared_ptr<SinusoidalTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    AbstractTimeDependence::setTimeDependence("td2", td2);
    AbstractTimeDependence::setTimeDependence("td3", td2);
    std::string name1 = AbstractTimeDependence::getName(td1);
    EXPECT_EQ(name1, "td1");
    std::string name2 = AbstractTimeDependence::getName(td2);
    EXPECT_TRUE(name2 == "td2" || name2 == "td3");
    EXPECT_THROW(AbstractTimeDependence::getName(td3),
                 GeneralClassicException);

}
