/*
 *  Copyright (c) 2014, Chris Rogers
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
#include "Algorithms/PolynomialTimeDependence.h"

#include "opal_test_utilities/SilenceTest.h"

TEST(PolynomialTimeDependenceTest, PolynomialTimeDependenceTest) {
    OpalTestUtilities::SilenceTest silencer;

    // Check empty polynomial coefficients always returns 0.
    std::vector<double> test;
    PolynomialTimeDependence time_dependence_1(test);
    EXPECT_DOUBLE_EQ(time_dependence_1.getValue(0.1), 0.);

    // Check constant term produces constant
    test.push_back(1.);
    PolynomialTimeDependence time_dependence_2(test);
    EXPECT_DOUBLE_EQ(time_dependence_2.getValue(0.1), 1.);

    // Check cubic terms
    test.push_back(2.);
    test.push_back(3.);
    PolynomialTimeDependence time_dependence_3(test);
    EXPECT_DOUBLE_EQ(time_dependence_3.getValue(0.1), 1.23);

    // Check clone produces same result
    PolynomialTimeDependence* time_dependence_clone = time_dependence_3.clone();
    EXPECT_DOUBLE_EQ(time_dependence_clone->getValue(0.1), 1.23);
    delete time_dependence_clone;
}

TEST(PolynomialTimeDependenceTest, TDMapTest) {
    OpalTestUtilities::SilenceTest silencer;

    // throw on empty value
    EXPECT_THROW(AbstractTimeDependence::getTimeDependence("name"),
                 GeneralClassicException);
    std::vector<double> test;

    // set/get time dependence
    PolynomialTimeDependence time_dep(test);
    std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td2", td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td2"), td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    // set time dependence overwriting existing time dependence
    // should overwrite, without memory leak
    std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td3);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td3);
}

TEST(PolynomialTimeDependenceTest, TDMapNameLookupTest) {
    OpalTestUtilities::SilenceTest silencer;

    EXPECT_THROW(AbstractTimeDependence::getName(NULL),
                 GeneralClassicException);
    PolynomialTimeDependence time_dep(std::vector<double>(1, 1));
    std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
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