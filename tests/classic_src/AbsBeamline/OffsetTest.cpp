//
// Unit tests for class Offset
//
// Copyright (c) 2014, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "opal_test_utilities/SilenceTest.h"

#include "AbsBeamline/Offset.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

#include "gtest/gtest.h"

#include <cmath>
#include <vector>

class SRotatedGeometry;

// see also RingSectionTest

TEST(OffsetTest, TestConstructDestruct) {
    OpalTestUtilities::SilenceTest silencer;

    // Test default constructors (copy constructor/operator below)
    Offset off1("Name");
    try {
      EXPECT_EQ(off1.getName(), "Name");
      EXPECT_TRUE(off1.getGeometry().getTotalTransform().isIdentity());

      Offset off2;
      EXPECT_EQ(off2.getName(), "");
      EXPECT_TRUE(off1.getGeometry().getTotalTransform().isIdentity());
    } catch(GeneralClassicException& exc) {
      std::cerr << exc.what() << std::endl;
    }
}

TEST(OffsetTest, TestGetSet) {
    OpalTestUtilities::SilenceTest silencer;

    // Check for typos in accessors/mutators
    Vector_t ref({1., 2., 3.});
    Vector_t test;
    Offset off;

    off.setEndDirection(ref);
    test = off.getEndDirection();
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(test(i), ref(i));
        ref(i) += 10;
    }

    off.setEndPosition(ref);
    test = off.getEndPosition();
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(test(i), ref(i));
        ref(i) += 10;
    }

    EXPECT_FALSE(off.getIsLocal());
    off.setIsLocal(true);
    EXPECT_TRUE(off.getIsLocal());
}

// check offset has desired input/output angles
void testOffset(Offset& off, double refRotIn, double refRotOut, double length,
                std::string msg) {
    Euclid3D transform = off.getGeometry().getTotalTransform();
    double rotIn = -std::atan2(transform.getVector()(0), transform.getVector()(2));
    double rotOut = -transform.getRotation().getAxis()(1);
    EXPECT_NEAR(refRotIn, rotIn, 1e-6) << msg+"\n";
    EXPECT_NEAR(refRotOut, rotOut, 1e-6) << msg+"\n";
    EXPECT_NEAR(length, off.getGeometry().getElementLength(), 1e-6) << msg+"\n";
}

// build an offset, call update geometry, check it was built
// supports updategeometry tests
void buildTestOffset(Vector_t endPos, Vector_t endDir,
                 double refRotIn, double refRotOut, double length,
                 std::string msg) {
    Offset off;
    off.setEndPosition(endPos);
    off.setEndDirection(endDir);
    off.setIsLocal(true);
    off.updateGeometry();
    testOffset(off, refRotIn, refRotOut, length, msg);
}

TEST(OffsetTest, TestUpdateIdentityTransforms) {
    OpalTestUtilities::SilenceTest silencer;

    // Check we don't throw up on an identity transform about origin
    // results are undefined, returns identity
    buildTestOffset(Vector_t({0., 0., 0.}), Vector_t({0., 1., 0.}),
                0., 0., 0., "identity transform about origin");

    // Check we don't throw up on a pure rotation
    buildTestOffset(Vector_t({0., 0., 0.}), Vector_t({-1., 0., 0.}),
                0., Physics::pi/2., 0., "Pure rotation thru pi/2");
}

TEST(OffsetTest, TestUpdateRotationsNotXY) {
    OpalTestUtilities::SilenceTest silencer;

    // Check we throw for a translation/rotation out of the midplane
    EXPECT_THROW(buildTestOffset(Vector_t({1., 3., 1.}), Vector_t({+1., 0., 0.}),
                             0., 0., 0., "not x-y rotation"),
                             GeneralClassicException);
    EXPECT_THROW(buildTestOffset(Vector_t({1., 3., 0.}), Vector_t({+1., 0., 1.}),
                             0., 0., 0., "not x-y rotation"),
                             GeneralClassicException);
}

TEST(OffsetTest, TestUpdateRotations) {
    OpalTestUtilities::SilenceTest silencer;

    // Check we get length right
    buildTestOffset(Vector_t({-2., 1., 0.}), Vector_t({+0., +1., 0.}),
                std::atan2(2., 1.), 0., std::sqrt(5.), "length");
    // theta_in is the rotation from start direction TO the displacement vector
    // where positive angle means rotation anticlockwise
    buildTestOffset(Vector_t({-1., 1., 0.}), Vector_t({0., 1., 0.}),
                +1.*Physics::pi/4., 0., std::sqrt(2.), "x-y rotation theta_in 1");
    buildTestOffset(Vector_t({-1., -1., 0.}), Vector_t({0., 1., 0.}),
                +3.*Physics::pi/4., 0., std::sqrt(2.), "x-y rotation theta_in 2");
    buildTestOffset(Vector_t({1., -1., 0.}), Vector_t({0., 1., 0.}),
                -3.*Physics::pi/4., 0., std::sqrt(2.), "x-y rotation theta_in 3");
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({0., 1., 0.}),
                -1.*Physics::pi/4., 0., std::sqrt(2.), "x-y rotation theta_in 4");
    // theta_out is the rotation from displacement vector TO end direction
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({0., 1., 0.}),
                -Physics::pi/4., 0., std::sqrt(2.),
                "x-y rotation theta_out 1");
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({-1., 1., 0.}),
                -Physics::pi/4., Physics::pi/4., std::sqrt(2.),
                "x-y rotation theta_out 2");
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({-1., 0., 0.}),
                -Physics::pi/4., 2.*Physics::pi/4., std::sqrt(2.),
                "x-y rotation theta_out 3");
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({1., -1., 0.}),
                -Physics::pi/4., -3.*Physics::pi/4., std::sqrt(2.),
                "x-y rotation theta_out 4");
    buildTestOffset(Vector_t({1., 1., 0.}), Vector_t({1., 1., 0.}),
                -Physics::pi/4., -1.*Physics::pi/4., std::sqrt(2.),
                "x-y rotation theta_out 5");
}

TEST(OffsetTest, TestCopy) {
    OpalTestUtilities::SilenceTest silencer;

    Offset off1("Name");
    off1.setEndPosition(Vector_t({1., 3., 0.}));
    off1.setEndDirection(Vector_t({1., 4., 0.}));
    off1.setIsLocal(true);
    off1.updateGeometry();
    Offset off2("Name", off1);
    // should be a deep copy
    EXPECT_NE(&off2.getGeometry(), &off1.getGeometry());
    // should be equal

    Offset off3;
    off3.setEndPosition(Vector_t({99., 99., 0.}));
    off3 = off2;
    EXPECT_NE(&off3.getGeometry(), &off2.getGeometry());
    EXPECT_TRUE(off1 == off3);
}

TEST(OffsetTest, TestRotateGetTheta) {
    OpalTestUtilities::SilenceTest silencer;

    // rotate and getTheta methods are (more or less) inverse
    // we check a few examples by hand
    Vector_t vecIn({1., 2., 0.});
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, vecIn), 0.);

    Vector_t rotOut = Offset::rotate(vecIn, Physics::pi/2.);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t({-2, 1, 0.})(i));
    }
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), Physics::pi/2.);

    vecIn = Vector_t({-1., -2., 0.});
    rotOut = Offset::rotate(vecIn, Physics::pi);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t({1, 2, 0.})(i));
    }
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), Physics::pi);

    vecIn = Vector_t({1., 2., 0.});
    rotOut = Offset::rotate(vecIn, 1.5*Physics::pi);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t({2, -1, 0.})(i));
    }
    // theta returns in domain -pi, pi
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), -0.5*Physics::pi);

    // we check many more examples automatically (using inverse property)
    for (double f = -Physics::pi; f < Physics::pi; f += Physics::pi/17.) {
        vecIn = Vector_t({1., 2., 0.});
        Vector_t rotOut1 = Offset::rotate(vecIn, f);
        Vector_t rotOut2 = Offset::rotate(vecIn, f+2.*Physics::pi);
        Vector_t rotOut3 = Offset::rotate(vecIn, f-2.*Physics::pi);
        EXPECT_NEAR(Offset::getTheta(vecIn, rotOut1), f, 1e-9) << "f=" << f;
        EXPECT_NEAR(Offset::getTheta(vecIn, rotOut2), f, 1e-9) << "f=" << f;
        EXPECT_NEAR(Offset::getTheta(vecIn, rotOut3), f, 1e-9) << "f=" << f;
    }
}

TEST(OffsetTest, TestBends) {
    OpalTestUtilities::SilenceTest silencer;
    double theta1 = Offset::float_tolerance*10.; // precision not great
    double theta2 = Offset::float_tolerance/1000.;
    Offset off = Offset::localCylindricalOffset("lco", theta1, 0., 3.);
    EXPECT_TRUE(off.bends());
    off = Offset::localCylindricalOffset("lco", 0., theta1, 3.);
    EXPECT_TRUE(off.bends());
    off = Offset::localCylindricalOffset("lco", theta1, -theta1, 3.);
    EXPECT_TRUE(off.bends());  // a chicane is considered to "bend"
    off = Offset::localCylindricalOffset("lco", theta2, theta2, 3.);
    EXPECT_FALSE(off.bends());
}

TEST(OffsetTest, TestLocalCylindricalOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/3.;
    Offset off1 = Offset::localCylindricalOffset("lco", theta, 0., 3.);
    EXPECT_EQ(off1.getName(), "lco");
    EXPECT_TRUE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t({-std::sin(theta), std::cos(theta), 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         Vector_t({-std::sin(theta), std::cos(theta), 0.})(i));

    Offset off2 = Offset::localCylindricalOffset("lco", 0., theta, 3.);
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off2.getEndPosition()(i), Vector_t({0., 3., 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off2.getEndDirection()(i),
                         Vector_t({-std::sin(theta), std::cos(theta), 0.})(i)) << i;

    Offset off3 = Offset::localCylindricalOffset("lco", theta, theta/3., 3.);
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off3.getEndPosition()(i),
                         3.*Vector_t({-std::sin(theta), std::cos(theta), 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off3.getEndDirection()(i),
                         Vector_t({-std::sin(4.*theta/3.), std::cos(4.*theta/3.), 0.})(i)) << i;
    //  testOffset(off1, theta, theta/3., 3., "");
}

TEST(OffsetTest, TestGlobalCylindricalOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double radius = 7.;
    double phi = Physics::pi/3.;
    double theta = Physics::pi/4.;
    Offset off1 = Offset::globalCylindricalOffset("gco", radius, phi, theta);
    EXPECT_EQ(off1.getName(), "gco");
    EXPECT_FALSE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         radius*Vector_t({std::cos(phi), std::sin(phi), 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         Vector_t({std::sin(theta+phi), std::cos(theta+phi), 0.})(i)) << i;
    testOffset(off1, 0., 0., 0., "");
}

TEST(OffsetTest, TestLocalCartesianOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/6.;
    Offset off1 = Offset::localCartesianOffset("lco",
                          3.*Vector_t({-std::sin(theta), std::cos(theta), 0.}),
                          10.*Vector_t({-std::sin(theta/3.), std::cos(theta/3.), 0.}));
    EXPECT_EQ(off1.getName(), "lco");
    EXPECT_TRUE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t({-std::sin(theta), std::cos(theta), 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         10.*Vector_t({-std::sin(theta/3.), std::cos(theta/3.), 0.})(i));
    testOffset(off1, theta, theta/3., 3., "");
}

TEST(OffsetTest, TestGlobalCartesianOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/3.;
    Offset off1 = Offset::globalCartesianOffset("gco",
                          3.*Vector_t({std::cos(theta), std::sin(theta), 0.}),
                          10.*Vector_t({std::cos(theta/3.), std::sin(theta/3.), 0.}));
    EXPECT_EQ(off1.getName(), "gco");
    EXPECT_FALSE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t({std::cos(theta), std::sin(theta), 0.})(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         10.*Vector_t({std::cos(theta/3.), std::sin(theta/3.), 0.})(i));
    testOffset(off1, 0., 0., 0., "");
}
