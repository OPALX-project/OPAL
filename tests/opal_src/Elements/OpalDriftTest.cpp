//
// Unit tests for OpalDrift element definition
//
// Copyright (c) 2026
// All rights reserved.
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
#include "gtest/gtest.h"

#include "Attributes/Attributes.h"
#include "Elements/OpalDrift.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

std::unique_ptr<OpalDrift> makeDrift(const std::optional<std::string>& aperture) {
    auto drift = std::make_unique<OpalDrift>();
    Attributes::setReal(*drift->findAttribute("L"), 2.0);
    if (aperture.has_value()) {
        Attributes::setString(*drift->findAttribute("APERTURE"), aperture.value());
    }
    drift->update();
    return drift;
}

}  // namespace

TEST(OpalDriftTest, CircleDefaultsMatchDefaultApertureBehaviour) {
    auto noAperture = makeDrift(std::nullopt);
    auto circle = makeDrift("CIRCLE(1)");
    auto conicCircle = makeDrift("CIRCLE(1,1)");

    ElementBase* noApertureElement = noAperture->getElement();
    ElementBase* circleElement = circle->getElement();
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(noApertureElement, nullptr);
    ASSERT_NE(circleElement, nullptr);
    ASSERT_NE(conicCircleElement, nullptr);

    const std::vector<Vector_t> probes = {
        Vector_t({0.00, 0.00, 0.20}),
        Vector_t({0.20, 0.10, 1.00}),
        Vector_t({0.49, 0.00, 0.20}),
        Vector_t({0.49, 0.00, 1.00}),
        Vector_t({0.49, 0.00, 1.80}),
        Vector_t({0.50, 0.00, 1.00}),
        Vector_t({0.00, 0.50, 1.00}),
        Vector_t({0.36, 0.36, 1.00}),
        Vector_t({0.40, 0.40, 1.00}),
        Vector_t({0.00, 0.00, -0.10}),
        Vector_t({0.00, 0.00, 2.00})
    };

    for (const Vector_t& r: probes) {
        SCOPED_TRACE(::testing::Message()
                     << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        EXPECT_EQ(noApertureElement->isInside(r), circleElement->isInside(r));
        EXPECT_EQ(noApertureElement->isInside(r), conicCircleElement->isInside(r));
    }
}

TEST(OpalDriftTest, SquareAndRectangleEquivalentBehaviour) {
    auto rectangle = makeDrift("RECTANGLE(1,1)");
    auto conicRectangle = makeDrift("RECTANGLE(1,1,1)");
    auto square = makeDrift("SQUARE(1)");
    auto conicSquare = makeDrift("SQUARE(1,1)");

    ElementBase* rectangleElement = rectangle->getElement();
    ElementBase* conicRectangleElement = conicRectangle->getElement();
    ElementBase* squareElement = square->getElement();
    ElementBase* conicSquareElement = conicSquare->getElement();

    ASSERT_NE(rectangleElement, nullptr);
    ASSERT_NE(conicRectangleElement, nullptr);
    ASSERT_NE(squareElement, nullptr);
    ASSERT_NE(conicSquareElement, nullptr);

    const std::vector<Vector_t> probes = {
        Vector_t({0.00, 0.00, 0.20}),
        Vector_t({0.49, 0.49, 1.00}),
        Vector_t({0.49, 0.10, 1.80}),
        Vector_t({0.10, 0.49, 0.20}),
        Vector_t({0.50, 0.10, 1.00}),
        Vector_t({0.10, 0.50, 1.00}),
        Vector_t({0.60, 0.40, 1.00}),
        Vector_t({0.40, 0.60, 1.00}),
        Vector_t({0.00, 0.00, -0.10}),
        Vector_t({0.00, 0.00, 2.00})
    };

    for (const Vector_t& r: probes) {
        SCOPED_TRACE(::testing::Message()
                     << "Probe point (" << r[0] << ", " << r[1] << ", " << r[2] << ")");
        const bool expected = rectangleElement->isInside(r);
        EXPECT_EQ(expected, conicRectangleElement->isInside(r));
        EXPECT_EQ(expected, squareElement->isInside(r));
        EXPECT_EQ(expected, conicSquareElement->isInside(r));
    }
}

TEST(OpalDriftTest, ConicCircleOpeningBehaviour) {
    auto conicCircle = makeDrift("CIRCLE(1,2)");
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(conicCircleElement, nullptr);

    const Vector_t centerline({0.00, 0.00, 1.00});
    EXPECT_TRUE(conicCircleElement->isInside(centerline));

    // Radius starts at 0.5 m and grows linearly to 1.0 m at the exit.
    const Vector_t startProbe({0.75, 0.00, 0.10});
    const Vector_t middleProbe({0.75, 0.00, 1.00});
    const Vector_t endProbe({0.75, 0.00, 1.90});

    EXPECT_FALSE(conicCircleElement->isInside(startProbe));
    EXPECT_FALSE(conicCircleElement->isInside(middleProbe));
    EXPECT_TRUE(conicCircleElement->isInside(endProbe));
}

TEST(OpalDriftTest, ConicCircleClosingBehaviour) {
    auto conicCircle = makeDrift("CIRCLE(1,0.5)");
    ElementBase* conicCircleElement = conicCircle->getElement();

    ASSERT_NE(conicCircleElement, nullptr);

    // Radius starts at 0.5 m and shrinks linearly to 0.25 m at the exit.
    const Vector_t startProbe({0.40, 0.00, 0.10});
    const Vector_t middleProbe({0.40, 0.00, 1.00});
    const Vector_t endProbe({0.40, 0.00, 1.90});

    EXPECT_TRUE(conicCircleElement->isInside(startProbe));
    EXPECT_FALSE(conicCircleElement->isInside(middleProbe));
    EXPECT_FALSE(conicCircleElement->isInside(endProbe));
}

TEST(OpalDriftTest, CircleConstantAlongLengthAndLongitudinalBounds) {
    auto circle = makeDrift("CIRCLE(1)");
    ElementBase* circleElement = circle->getElement();

    ASSERT_NE(circleElement, nullptr);

    // Same transverse position should give the same result at multiple z inside the element.
    const std::vector<Vector_t> insideTransverseInsideLength = {
        Vector_t({0.49, 0.00, 0.10}),
        Vector_t({0.49, 0.00, 1.00}),
        Vector_t({0.49, 0.00, 1.90})
    };
    for (const Vector_t& r: insideTransverseInsideLength) {
        EXPECT_TRUE(circleElement->isInside(r));
    }

    const std::vector<Vector_t> outsideTransverseInsideLength = {
        Vector_t({0.51, 0.00, 0.10}),
        Vector_t({0.51, 0.00, 1.00}),
        Vector_t({0.51, 0.00, 1.90})
    };
    for (const Vector_t& r: outsideTransverseInsideLength) {
        EXPECT_FALSE(circleElement->isInside(r));
    }

    // Before and after the element in z should always be outside, even on axis.
    EXPECT_FALSE(circleElement->isInside(Vector_t({0.00, 0.00, -0.10})));
    EXPECT_FALSE(circleElement->isInside(Vector_t({0.00, 0.00, 2.00})));
}
