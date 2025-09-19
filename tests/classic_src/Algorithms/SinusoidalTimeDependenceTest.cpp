//
// Tests for SinusoidalTimeDependence
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

TEST(SinusoidalTimeDependenceTest, Integral) {
    OpalTestUtilities::SilenceTest silencer;

    // Check empty coefficients always returns 0.
    SinusoidalTimeDependence time_dependence_1({}, {}, {}, {});
    EXPECT_DOUBLE_EQ(time_dependence_1.getIntegral(0.1), 0.0);

    // Check a sine wave value
    SinusoidalTimeDependence time_dependence_2({8.0}, {}, {}, {});
    EXPECT_DOUBLE_EQ(time_dependence_2.getIntegral(0.1), 0.013746670117215259);

    // Check the amplitude
    SinusoidalTimeDependence time_dependence_3({8.0}, {}, {2.0}, {});
    EXPECT_DOUBLE_EQ(time_dependence_3.getIntegral(0.1), 0.013746670117215259 * 2);

    // Check the phase offset
    SinusoidalTimeDependence time_dependence_4({8.0}, {0.1}, {2.0} , {});
    EXPECT_DOUBLE_EQ(time_dependence_4.getIntegral(0.1), 0.02357815814417235);

    // Check the DC offset
    SinusoidalTimeDependence time_dependence_5({8.0}, {0.1}, {2.0}, {-1.0});
    EXPECT_DOUBLE_EQ(time_dependence_5.getIntegral(0.1), 0.02357815814417235 - 0.1);
}
