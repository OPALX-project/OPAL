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

#include <vector>

#include "gtest/gtest.h"

#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "AbsBeamline/Offset.h"

#include "opal_test_utilities/SilenceTest.h"

class SRotatedGeometry;

const double units = 1e3;

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
    Vector_t ref(1., 2., 3.);
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
    double rotIn = -atan2(transform.getVector()(0), transform.getVector()(2));
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
    buildTestOffset(Vector_t(0., 0., 0.), Vector_t(1., 0., 0.),
                0., 0., 0., "identity transform about origin");

    // Check we don't throw up on a pure rotation
    buildTestOffset(Vector_t(0., 0., 0.), Vector_t(0., 1., 0.),
                0., Physics::pi/2., 0., "Pure rotation thru pi/2");
}

TEST(OffsetTest, TestUpdateRotationsNotXY) {
    OpalTestUtilities::SilenceTest silencer;

    // Check we throw for a translation/rotation out of the midplane
    EXPECT_THROW(buildTestOffset(Vector_t(1., 3., 1.), Vector_t(+1., 0., 0.),
                             0., 0., 0., "not x-y rotation"),
                             GeneralClassicException);
    EXPECT_THROW(buildTestOffset(Vector_t(1., 3., 0.), Vector_t(+1., 0., 1.),
                             0., 0., 0., "not x-y rotation"),
                             GeneralClassicException);
}

TEST(OffsetTest, TestUpdateRotations) {
    OpalTestUtilities::SilenceTest silencer;

    // Check we get length right
    buildTestOffset(Vector_t(1., 2., 0.), Vector_t(+1., 0., 0.),
                atan2(2., 1.), 0., sqrt(5.), "length");
    // theta_in is the rotation from start direction TO the displacement vector
    // displacement vector is 1., 1., 0.
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(1., 0., 0.),
                +1.*Physics::pi/4., 0., sqrt(2.), "x-y rotation theta_in 1");
    buildTestOffset(Vector_t(-1., 1., 0.), Vector_t(1., 0., 0.),
                +3.*Physics::pi/4., 0., sqrt(2.), "x-y rotation theta_in 2");
    buildTestOffset(Vector_t(-1., -1., 0.), Vector_t(1., 0., 0.),
                -3.*Physics::pi/4., 0., sqrt(2.), "x-y rotation theta_in 3");
    buildTestOffset(Vector_t(1., -1., 0.), Vector_t(1., 0., 0.),
                -1.*Physics::pi/4., 0., sqrt(2.), "x-y rotation theta_in 4");
    // theta_out is the rotation from displacement vector TO end direction
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(1., 0., 0.),
                Physics::pi/4., 0., sqrt(2.),
                "x-y rotation theta_out 1");
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(1., 1., 0.),
                Physics::pi/4., Physics::pi/4., sqrt(2.),
                "x-y rotation theta_out 2");
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(0., 1., 0.),
                Physics::pi/4., 2.*Physics::pi/4., sqrt(2.),
                "x-y rotation theta_out 3");
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(-1., -1., 0.),
                Physics::pi/4., -3.*Physics::pi/4., sqrt(2.),
                "x-y rotation theta_out 4");
    buildTestOffset(Vector_t(1., 1., 0.), Vector_t(1., -1., 0.),
                Physics::pi/4., -1.*Physics::pi/4., sqrt(2.),
                "x-y rotation theta_out 5");
}

TEST(OffsetTest, TestCopy) {
    OpalTestUtilities::SilenceTest silencer;

    Offset off1("Name");
    off1.setEndPosition(Vector_t(1., 3., 0.));
    off1.setEndDirection(Vector_t(1., 4., 0.));
    off1.setIsLocal(true);
    off1.updateGeometry();
    Offset off2("Name", off1);
    // should be a deep copy
    EXPECT_NE(&off2.getGeometry(), &off1.getGeometry());
    // should be equal

    Offset off3;
    off3.setEndPosition(Vector_t(99., 99., 0.));
    off3 = off2;
    EXPECT_NE(&off3.getGeometry(), &off2.getGeometry());
    EXPECT_TRUE(off1 == off3);
}

TEST(OffsetTest, TestRotateGetTheta) {
    OpalTestUtilities::SilenceTest silencer;

    // rotate and getTheta methods are (more or less) inverse
    // we check a few examples by hand
    Vector_t vecIn(1., 2., 0.);
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, vecIn), 0.);

    Vector_t rotOut = Offset::rotate(vecIn, Physics::pi/2.);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t(-2, 1, 0.)(i));
    }
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), Physics::pi/2.);

    vecIn = Vector_t(-1., -2., 0.);
    rotOut = Offset::rotate(vecIn, Physics::pi);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t(1, 2, 0.)(i));
    }
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), Physics::pi);

    vecIn = Vector_t(1., 2., 0.);
    rotOut = Offset::rotate(vecIn, 1.5*Physics::pi);
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(rotOut(i), Vector_t(2, -1, 0.)(i));
    }
    // theta returns in domain -pi, pi
    EXPECT_DOUBLE_EQ(Offset::getTheta(vecIn, rotOut), -0.5*Physics::pi);

    // we check many more examples automatically (using inverse property)
    for (double f = -Physics::pi; f < Physics::pi; f += Physics::pi/17.) {
        vecIn = Vector_t(1., 2., 0.);
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
    Offset off = Offset::localCylindricalOffset("lco", theta1, 0., 3./units);
    EXPECT_TRUE(off.bends());
    off = Offset::localCylindricalOffset("lco", 0., theta1, 3./units);
    EXPECT_TRUE(off.bends());
    off = Offset::localCylindricalOffset("lco", theta1, -theta1, 3./units);
    EXPECT_TRUE(off.bends());  // a chicane is considered to "bend"
    off = Offset::localCylindricalOffset("lco", theta2, theta2, 3./units);
    EXPECT_FALSE(off.bends());
}

TEST(OffsetTest, TestLocalCylindricalOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/3.;
    Offset off1 = Offset::localCylindricalOffset("lco", theta, 0., 3./units);
    EXPECT_EQ(off1.getName(), "lco");
    EXPECT_TRUE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t(cos(theta), sin(theta), 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         Vector_t(cos(theta), sin(theta), 0.)(i));

    Offset off2 = Offset::localCylindricalOffset("lco", 0., theta, 3./units);
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off2.getEndPosition()(i), Vector_t(3., 0., 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off2.getEndDirection()(i),
                         Vector_t(cos(theta), sin(theta), 0.)(i)) << i;

    Offset off3 = Offset::localCylindricalOffset("lco", theta, theta/3., 3./units);
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off3.getEndPosition()(i),
                         3.*Vector_t(cos(theta), sin(theta), 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off3.getEndDirection()(i),
                         Vector_t(cos(4.*theta/3.), sin(4.*theta/3.), 0.)(i)) << i;
    //  testOffset(off1, theta, theta/3., 3., "");
}

TEST(OffsetTest, TestGlobalCylindricalOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double radius = 7.;
    double phi = Physics::pi/3.;
    double theta = Physics::pi/4.;
    Offset off1 = Offset::globalCylindricalOffset("gco", radius/units, phi, theta);
    EXPECT_EQ(off1.getName(), "gco");
    EXPECT_FALSE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         radius*Vector_t(cos(phi), sin(phi), 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         Vector_t(sin(theta+phi), cos(theta+phi), 0.)(i)) << i;
    testOffset(off1, 0., 0., 0., "");
}

TEST(OffsetTest, TestLocalCartesianOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/6.;
    Offset off1 = Offset::localCartesianOffset("lco",
                          3./units*Vector_t(cos(theta), sin(theta), 0.),
                          10.*Vector_t(cos(theta/3.), sin(theta/3.), 0.));
    EXPECT_EQ(off1.getName(), "lco");
    EXPECT_TRUE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t(cos(theta), sin(theta), 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         10.*Vector_t(cos(theta/3.), sin(theta/3.), 0.)(i));
    testOffset(off1, theta, theta/3., 3., "");
}

TEST(OffsetTest, TestGlobalCartesianOffset) {
    OpalTestUtilities::SilenceTest silencer;

    double theta = Physics::pi/3.;
    Offset off1 = Offset::globalCartesianOffset("gco",
                          3./units*Vector_t(cos(theta), sin(theta), 0.),
                          10.*Vector_t(cos(theta/3.), sin(theta/3.), 0.));
    EXPECT_EQ(off1.getName(), "gco");
    EXPECT_FALSE(off1.getIsLocal());
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndPosition()(i),
                         3.*Vector_t(cos(theta), sin(theta), 0.)(i)) << i;
    for (int i = 0; i<3; ++i)
        EXPECT_DOUBLE_EQ(off1.getEndDirection()(i),
                         10.*Vector_t(cos(theta/3.), sin(theta/3.), 0.)(i));
    testOffset(off1, 0., 0., 0., "");
}