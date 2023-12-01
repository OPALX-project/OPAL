// Copyright (c) 2023, Chris Rogers
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


#include "gtest/gtest.h"

#include "Attributes/Attributes.h"
#include "AbsBeamline/Offset.h"
#include "Elements/OpalOutputPlane.h"
#include "AbsBeamline/OutputPlane.h"

TEST(OpalOutputPlaneTest, TestCtor) {
    OpalOutputPlane plane1;
    EXPECT_EQ((&plane1)->getOpalName(), "OUTPUTPLANE");
}


void setReal(OpalOutputPlane& plane, std::string attName, double value) {
    Attribute* att = plane.findAttribute(attName);
    ASSERT_NE(att, nullptr) << "Failed with attribute " << attName;
    Attributes::setReal(*att, value);
}

TEST(OpalOutputPlaneTest, TestProbeMode) {
    OpalOutputPlane plane;
    Attribute* att = plane.findAttribute("PLACEMENT_STYLE");
    Attributes::setPredefinedString(*att, "PROBE");
    setReal(plane, "XSTART", 1.0);
    setReal(plane, "XEND", 3.0);
    setReal(plane, "YSTART", -1.0);
    setReal(plane, "YEND", -3.0);
    plane.update();
    OutputPlane *output = dynamic_cast<OutputPlane *>(plane.getElement());
    double roottwo = std::sqrt(2);
    EXPECT_NEAR(output->getCentre()[0], 2.0, 1e-12);
    EXPECT_NEAR(output->getCentre()[1], -2.0, 1e-12);
    EXPECT_NEAR(output->getCentre()[2], 0.0, 1e-12);
    EXPECT_NEAR(output->getHorizontalExtent(), roottwo, 1e-12); // that is half width
    EXPECT_NEAR(output->getNormal()[0], 1/roottwo, 1e-12);
    EXPECT_NEAR(output->getNormal()[1], 1/roottwo, 1e-12);

    // swap and check normal is reversed
    setReal(plane, "XSTART", 3.0);
    setReal(plane, "XEND", 1.0);
    setReal(plane, "YSTART", -3.0);
    setReal(plane, "YEND", -1.0);
    plane.update();
    EXPECT_NEAR(output->getNormal()[0], -1/roottwo, 1e-12);
    EXPECT_NEAR(output->getNormal()[1], -1/roottwo, 1e-12);
}

TEST(OpalOutputPlaneTest, TestPlacement) {
    OpalOutputPlane plane;
    Attribute* att0 = plane.findAttribute("PLACEMENT_STYLE");
    Attributes::setPredefinedString(*att0, "CENTRE_NORMAL");
    Attribute* att1 = plane.findAttribute("CENTRE");
    Attributes::setRealArray(*att1, {1.0, 2.0, 3.0});
    Attribute* att2 = plane.findAttribute("NORMAL");
    Attributes::setRealArray(*att2, {-4.0, 4.0, 2.0});
    setReal(plane, "WIDTH", 4.0);
    setReal(plane, "HEIGHT", 5.0);
    setReal(plane, "RADIUS", 6.0);
    plane.update();
    OutputPlane *output = dynamic_cast<OutputPlane *>(plane.getElement());
    EXPECT_NEAR(output->getCentre()[0], 1.0, 1e-12);
    EXPECT_NEAR(output->getCentre()[1], 2.0, 1e-12);
    EXPECT_NEAR(output->getCentre()[2], 3.0, 1e-12);
    EXPECT_NEAR(output->getNormal()[0], -2/3.0, 1e-12);
    EXPECT_NEAR(output->getNormal()[1], 2/3.0, 1e-12);
    EXPECT_NEAR(output->getNormal()[2], 1.0/3.0, 1e-12);
    EXPECT_NEAR(output->getHorizontalExtent(), 2.0, 1e-12); // that is half width
    EXPECT_NEAR(output->getVerticalExtent(), 2.5, 1e-12); // that is half height
    EXPECT_NEAR(output->getRadialExtent(), 6.0, 1e-12); // that is half height
}
