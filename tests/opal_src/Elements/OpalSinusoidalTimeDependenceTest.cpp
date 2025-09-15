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
#include <Algorithms/AbstractTimeDependence.h>
#include <Attributes/Attributes.h>

#include <sstream>

#include "gtest/gtest.h"
#include "Elements/OpalSinusoidalTimeDependence.h"
#include "opal_test_utilities/SilenceTest.h"

class BeamlineVisitor;
// some comment
class TestElement final : public ElementBase {
public:

    TestElement() : base(nullptr) {}
    ElementType getType() const override {
        return ElementType::ANY;}
    BGeometryBase  &getGeometry() override {
        return *base;
    }

    const BGeometryBase  &getGeometry() const override {
        return *base;
    }

    ElementBase* clone() const override {return nullptr;}

    void accept(BeamlineVisitor& /*visitor*/) const override {}

private:
    BGeometryBase* base;
    std::string type{};
};

TEST(OpalSinusoidalTimeDependenceTest, ConstructorTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalSinusoidalTimeDependence dep;
    const OpalSinusoidalTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(OpalSinusoidalTimeDependenceTest, PrintTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalSinusoidalTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "SINUSOIDAL_TIME_DEPENDENCE;\n");
}

TEST(OpalSinusoidalTimeDependenceTest, UpdateTest) {
    OpalTestUtilities::SilenceTest silencer;
    OpalSinusoidalTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setRealArray(dependence.itsAttr[OpalSinusoidalTimeDependence::FREQUENCIES],
        {8.0});
    Attributes::setRealArray(dependence.itsAttr[OpalSinusoidalTimeDependence::AMPLITUDES],
        {3.0});
    Attributes::setRealArray(dependence.itsAttr[OpalSinusoidalTimeDependence::PHASE_OFFSETS],
        {0.1});
    Attributes::setRealArray(dependence.itsAttr[OpalSinusoidalTimeDependence::DC_OFFSETS],
        {-1.0});
    dependence.update();
    const auto dep = AbstractTimeDependence::getTimeDependence("SCALE1");
    const auto value = dep->getValue(0.1);
    const auto shouldBe = 3.0 / 2.0 * std::sin(2.0 * M_PI * 8.0 * 0.1 + 0.1) - 1.0;
    EXPECT_DOUBLE_EQ(value, shouldBe);
}
