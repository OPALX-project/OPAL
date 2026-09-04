//
// Unit tests for OpalScalingFFAMagnet element definition
//
// Copyright (c) 2022, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "opal_test_utilities/SilenceTest.h"

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "gtest/gtest.h"

#include <cmath>
#include <sstream>

void setReal(OpalScalingFFAMagnet& mag, std::string attName, double value) {
    Attribute* att = mag.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

TEST(OpalScalingFFAMagnetTest, TestConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    OpalScalingFFAMagnet mag1;
    EXPECT_EQ(mag1.getOpalName(), "SCALINGFFAMAGNET");
}

TEST(OpalScalingFFAMagnetTest, TestUpdate) {
    OpalScalingFFAMagnet opalMag1;
    setReal(opalMag1, "B0", 1.0);
    setReal(opalMag1, "R0", 12.5);
    setReal(opalMag1, "FIELD_INDEX", 2.0);
    setReal(opalMag1, "TAN_DELTA", 0.7);
    setReal(opalMag1, "MAX_Y_POWER", 3.0);
    setReal(opalMag1, "END_LENGTH", 0.8);
    setReal(opalMag1, "CENTRE_LENGTH", 5.0);
    setReal(opalMag1, "RADIAL_NEG_EXTENT", 0.55);
    setReal(opalMag1, "RADIAL_POS_EXTENT", 0.65);
    setReal(opalMag1, "HEIGHT", 6.0);
    opalMag1.update();
    ElementBase* element = opalMag1.getElement();
    ScalingFFAMagnet* mag1 = dynamic_cast<ScalingFFAMagnet*>(element);
    mag1->setupEndField();
    ASSERT_NE(mag1, nullptr);
    EXPECT_NEAR(mag1->getDipoleConstant(), 1.0*Units::T2kG, 1e-9);
    EXPECT_NEAR(mag1->getR0(), 12.5, 1e-9);
    EXPECT_NEAR(mag1->getFieldIndex(), 2.0, 1e-9);
    EXPECT_NEAR(mag1->getTanDelta(), 0.7, 1e-9);
    EXPECT_EQ(mag1->getMaxOrder(), 3u);
    EXPECT_NEAR(mag1->getRMin(), 11.95, 1e-9);
    EXPECT_NEAR(mag1->getRMax(), 13.15, 1e-9);
    EXPECT_NEAR(mag1->getVerticalExtent(), 3.0, 1e-9);

    endfieldmodel::EndFieldModel* model = mag1->getEndField();
    endfieldmodel::Tanh* tanh = dynamic_cast<endfieldmodel::Tanh*>(model);
    ASSERT_NE(tanh, nullptr);
    EXPECT_NEAR(tanh->getX0(), 5.0/12.5/2, 1e-9);
    EXPECT_NEAR(tanh->getLambda(), 0.8/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiStart(), 0.5*(4*0.8+5.0)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiEnd(), (4*0.8+5.0)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getAzimuthalExtent(), (5*0.8+5.0/2.0)/12.5, 1e-9);

    setReal(opalMag1, "MAGNET_START", 7.0);
    setReal(opalMag1, "MAGNET_END", 8.0);
    setReal(opalMag1, "AZIMUTHAL_EXTENT", 9.0);
    opalMag1.update();
    mag1->setupEndField();
    EXPECT_NEAR(mag1->getPhiStart(), (7.0+5.0/2)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiEnd(), 8.0/12.5, 1e-9);
    EXPECT_NEAR(mag1->getAzimuthalExtent(), 9.0/12.5, 1e-9);
}

// r in mm, phi in rad
Vector_t cartesianCoord(double r0, double phi) {
    //phi -= Physics::pi/2.0;
    return Vector_t({r0*std::cos(phi)-r0, 0.0, r0*std::sin(phi)});
}

TEST(OpalScalingFFAMagnetTest, TestFieldCheck) {
    OpalScalingFFAMagnet opalMag1;
    double r0 = 4.0;
    double b0 = 1.1;
    setReal(opalMag1, "B0", b0);
    setReal(opalMag1, "R0", r0);
    setReal(opalMag1, "FIELD_INDEX", 2.0);
    setReal(opalMag1, "MAX_Y_POWER", 3.0);
    setReal(opalMag1, "END_LENGTH", r0*Physics::pi/64.0);
    setReal(opalMag1, "CENTRE_LENGTH", r0*16*Physics::pi/64.0); // 1/8 of a circle
    setReal(opalMag1, "RADIAL_NEG_EXTENT", 2.0);
    setReal(opalMag1, "RADIAL_POS_EXTENT", 2.0);
    setReal(opalMag1, "HEIGHT", 6.0);
    setReal(opalMag1, "MAGNET_START", r0*12*Physics::pi/64.0);
    setReal(opalMag1, "MAGNET_END", r0*32*Physics::pi/64.0);
    setReal(opalMag1, "AZIMUTHAL_EXTENT", r0*16.0*Physics::pi/64.0);
    opalMag1.update();
    ElementBase* element = opalMag1.getElement();
    ScalingFFAMagnet* mag1 = dynamic_cast<ScalingFFAMagnet*>(element);
    mag1->setupEndField();
    ASSERT_NE(mag1, nullptr);
    // I just want to check that the bounding boxes/etc make sense
    double b0kG = b0*Units::T2kG;
    std::vector<double> position = {3.99, 4.01, 12, 20, 28, 35.99, 36.01};
    std::vector<bool> oob = {1, 0, 0, 0, 0, 0, 1};
    std::vector<double> field = {0.0, 0.0, b0kG/2, b0kG, b0kG/2, 0.0, 0.0};
    for (size_t i = 0; i < position.size(); ++i) {
        double phi = position[i]*Physics::pi/64.0;
        Vector_t B, P, E;
        double t = 0.;
        // magnet start plus half centre length should have maximum field
        Vector_t rMiddle = cartesianCoord(r0, phi);
        bool outOfBounds = mag1->apply(rMiddle, P, t, E, B);
        EXPECT_EQ(outOfBounds, oob[i]);
        EXPECT_NEAR(B[1], field[i], 1e-4) << "failed for phi "
                                        << position[i] << "*PI/64 " << std::endl;
    }

    Euclid3D delta = mag1->getGeometry().getTotalTransform();
    Vector3D vec = delta.getVector();
    Vector3D rot = delta.getRotation().getAxis();
    EXPECT_NEAR(vec(0), -r0*std::sin(mag1->getPhiEnd()), 1e-6);
    EXPECT_NEAR(vec(1), 0., 1e-6);
    EXPECT_NEAR(vec(2), r0*(1-std::cos(mag1->getPhiEnd())), 1e-6);

    EXPECT_EQ(rot(0), 0.);
    EXPECT_EQ(rot(1), -mag1->getPhiEnd());
    EXPECT_EQ(rot(2), 0.);

    ////// REVERSE BEND ANGLE ///////////
    setReal(opalMag1, "R0", -r0);
    opalMag1.update();
    element = opalMag1.getElement();
    mag1 = dynamic_cast<ScalingFFAMagnet*>(element);
    mag1->setupEndField();

    Euclid3D delta2 = mag1->getGeometry().getTotalTransform();
    Vector3D vec2 = delta2.getVector();
    Vector3D rot2 = delta2.getRotation().getAxis();
    EXPECT_NEAR(vec2(0), r0*std::sin(mag1->getPhiEnd()), 1e-6);
    EXPECT_NEAR(vec2(1), 0., 1e-6);
    EXPECT_NEAR(vec2(2), r0*(1-std::cos(mag1->getPhiEnd())), 1e-6);

    EXPECT_EQ(rot2(0), 0.);
    EXPECT_EQ(rot2(1), mag1->getPhiEnd());
    EXPECT_EQ(rot2(2), 0.);
}
