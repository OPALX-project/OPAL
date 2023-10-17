//
// Unit tests for class RingSection
//
// Copyright (c) 2014 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "opal_src/Utilities/MockComponent.h"
#include "opal_test_utilities/SilenceTest.h"

#include "AbsBeamline/Offset.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/RingSection.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>

TEST(RingSectionTest, TestConstructDestruct) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    MockComponent* compNull = nullptr;
    Vector_t vec0(0, 0, 0);
    EXPECT_EQ(ors.getComponent(), compNull);
    EXPECT_EQ(ors.getStartPosition(), vec0);
    EXPECT_EQ(ors.getStartNormal(), vec0);
    EXPECT_EQ(ors.getEndPosition(), vec0);
    EXPECT_EQ(ors.getEndNormal(), vec0);
    EXPECT_EQ(ors.getComponentPosition(), vec0);
    EXPECT_EQ(ors.getComponentOrientation(), vec0);
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_EQ(ors.getVirtualBoundingBox().at(i)[j], 0.);

    RingSection ors_comp;
    MockComponent comp;
    ors_comp.setComponent(&comp);
    // and implicit destructors; should not double free comp;
}

TEST(RingSectionTest, TestIsOnOrPastStartPlane) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    ors.setStartPosition(Vector_t(0., 1., 0.));
    ors.setStartNormal(Vector_t(1., 0., 0.));
    Vector_t vec1(1e-9, 1.e-9, 0.);
    Vector_t vec2(-1e-9, 1.e-9, 0.);
    Vector_t vec3(1e-9, -1.e-9, 0.); // other side of the ring
    Vector_t vec4(-1e-9, -1.e-9, 0.); // other side of the ring
    Vector_t vec5(1e-9, 1.e9, 0.); // large radius
    Vector_t vec6(-1e-9, 1.e9, 0.); // large radius
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec1));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec2));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec3));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec4));
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec5));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec6));

    ors.setStartNormal(Vector_t(-1., 0., 0.)); // rotate 180 degrees
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec1));
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec2));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec3));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec4));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec5));
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec6));

    ors.setStartPosition(Vector_t(-1., -1., 0.));
    ors.setStartNormal(Vector_t(1., -0.5, 0.));

    Vector_t vec7(-1.1e-9, 1.e-9, 0.); // this side of the ring
    Vector_t vec8(-0.9e-9, 1.e-9, 0.); // other side of the ring
    Vector_t vec9(-0.5-1e-9, 0., 0.); // behind normal
    Vector_t vec10(-0.5+1e-9, 0., 0.); // in front of normal
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec7));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec8));
    EXPECT_FALSE(ors.isOnOrPastStartPlane(vec9));
    EXPECT_TRUE(ors.isOnOrPastStartPlane(vec10));
}

TEST(RingSectionTest, TestIsPastEndPlane) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    ors.setEndPosition(Vector_t(0., 1., 0.));
    ors.setEndNormal(Vector_t(1., 0., 0.));
    Vector_t vec1(1e-9, 1.e-9, 0.);
    Vector_t vec2(-1e-9, 1.e-9, 0.);
    Vector_t vec3(1e-9, -1.e-9, 0.); // other side of the ring
    Vector_t vec4(-1e-9, -1.e-9, 0.); // other side of the ring
    Vector_t vec5(1e-9, 1.e9, 0.); // large radius
    Vector_t vec6(-1e-9, 1.e9, 0.); // large radius
    EXPECT_TRUE(ors.isPastEndPlane(vec1));
    EXPECT_FALSE(ors.isPastEndPlane(vec2));
    EXPECT_FALSE(ors.isPastEndPlane(vec3));
    EXPECT_FALSE(ors.isPastEndPlane(vec4));
    EXPECT_TRUE(ors.isPastEndPlane(vec5));
    EXPECT_FALSE(ors.isPastEndPlane(vec6));

    ors.setEndNormal(Vector_t(-1., 0., 0.)); // rotate 180 degrees
    EXPECT_FALSE(ors.isPastEndPlane(vec1));
    EXPECT_TRUE(ors.isPastEndPlane(vec2));
    EXPECT_FALSE(ors.isPastEndPlane(vec3));
    EXPECT_FALSE(ors.isPastEndPlane(vec4));
    EXPECT_FALSE(ors.isPastEndPlane(vec5));
    EXPECT_TRUE(ors.isPastEndPlane(vec6));

    ors.setEndPosition(Vector_t(-1., -1., 0.));
    ors.setEndNormal(Vector_t(1., -0.5, 0.));

    Vector_t vec7(-1.1e-9, 1.e-9, 0.); // this side of the ring
    Vector_t vec8(-0.9e-9, 1.e-9, 0.); // other side of the ring
    Vector_t vec9(-0.5-1e-9, 0., 0.); // behind normal
    Vector_t vec10(-0.5+1e-9, 0., 0.); // in front of normal
    EXPECT_TRUE(ors.isPastEndPlane(vec7));
    EXPECT_FALSE(ors.isPastEndPlane(vec8));
    EXPECT_FALSE(ors.isPastEndPlane(vec9));
    EXPECT_TRUE(ors.isPastEndPlane(vec10));
}

TEST(RingSectionTest, TestGetFieldValue) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    MockComponent comp;
    ors.setComponent(&comp);
    Vector_t centre; //(-1.33, +1.66, 0.);
    for (double theta = -3.*Physics::pi; theta < 3.*Physics::pi; theta += Physics::pi/6.) {
        Vector_t orientation(0., 0., theta);
        ors.setComponentOrientation(orientation);
        ors.setComponentPosition(centre);
        double c = std::cos(orientation(2));
        double s = std::sin(orientation(2));
        // x y z are coordinates in local OPAL-CYCL coordinate system
        for (double x = 0.01; x < 1.; x += 0.1)
            for (double y = 0.01; y < 1.; y += 0.1)
                for (double z = -0.01; z > -1.; z -= 0.1) {
                    Vector_t offset(c*x+s*y, +s*x-c*y, z);
                    Vector_t pos = centre+offset;
                    Vector_t centroid, B, E;
                    double t = 0;
                    bool oobRef =  (x < 0. || x > 1. ||
                                    z < -1. || z > 0. ||
                                    y < 0. || y > 1.);
                    bool oobTest = ors.getFieldValue(pos, centroid, t, E, B);
                    EXPECT_EQ(oobTest, oobRef);
                    Vector_t bfield(c*x+s*y, s*x-c*y, z);
                    std::cout << "loc " << x << " " << y << " " << z 
                              << " glob " << pos 
                              << " bglob " << B 
                              << " bloc " << bfield 
                              << " oobTest " << oobTest 
                              << " oobRef " << oobRef << std::endl;
                    for (int l = 0; l < 3; ++l) {
                        EXPECT_NEAR(B(l), +bfield(l), 1e-6);
                        EXPECT_NEAR(E(l), -bfield(l), 1e-6);
                    }
                }
    }
}

namespace {
    bool sort_comparator(Vector_t v1, Vector_t v2) {
        if (std::abs(v1(0) - v2(0)) < 1e-6) {
            if (std::abs(v1(1) - v2(1)) < 1e-6) {
                return v1(2) > v2(2);
            }
            return v1(1) > v2(1);
        }
        return v1(0) > v2(0);
    }
}

TEST(RingSectionTest, TestGetVirtualBoundingBox) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    ors.setStartPosition(Vector_t(3, -1, 99));
    ors.setStartNormal(Vector_t(-4, -1, -1000));
    ors.setEndPosition(Vector_t(2, 1, 77));
    ors.setEndNormal(Vector_t(-1, 1, 655));
    std::vector<Vector_t> bb = ors.getVirtualBoundingBox();
    std::vector<Vector_t> bbRef;
    bbRef.push_back(Vector_t(0.99*std::sqrt(10)/(-std::sqrt(17))+3.,
                             0.99*std::sqrt(10)*4./(+std::sqrt(17))-1., 99.));
    bbRef.push_back(Vector_t(0.99*std::sqrt(10)/(+std::sqrt(17))+3.,
                             0.99*std::sqrt(10)*4./(-std::sqrt(17))-1., 99.));
    bbRef.push_back(Vector_t(0.99*std::sqrt(5)/(+std::sqrt(2))+2.,
                             0.99*std::sqrt(5)/(+std::sqrt(2))+1., 77.));
    bbRef.push_back(Vector_t(0.99*std::sqrt(5)/(-std::sqrt(2))+2.,
                             0.99*std::sqrt(5)/(-std::sqrt(2))+1., 77.));
    std::sort(bb.begin(), bb.end(), sort_comparator);
    std::sort(bbRef.begin(), bbRef.end(), sort_comparator);
    EXPECT_EQ(bb.size(), bbRef.size());
    for (size_t i = 0; i < bb.size(); ++i) {
        for (size_t j = 0; j < 3; ++j)
            EXPECT_NEAR(bb[i](j), bbRef[i](j), 1e-6);
    }
}

namespace {
    RingSection buildORS(double r, double phi1, double phi2) {
        RingSection ors;
        ors.setStartPosition(Vector_t(std::sin(phi1)*r, std::cos(phi1)*r, 0.));
        ors.setStartNormal(Vector_t(std::cos(phi1), -std::sin(phi1), 0.));
        ors.setEndPosition(Vector_t(std::sin(phi2)*r, std::cos(phi2)*r, 0.));
        ors.setEndNormal(Vector_t(std::cos(phi2), -std::sin(phi2), 0.));
        return ors;
    }
}

TEST(RingSectionTest, TestDoesOverlap) {
    OpalTestUtilities::SilenceTest silencer;

    double f1 = 1.0*Physics::pi/6.;
    double f2 = 0.5*Physics::pi/6.;
    double f3 = -0.5*Physics::pi/6.;
    double f4 = -1.0*Physics::pi/6.;
    double r = 3.;
    RingSection ors1 = buildORS(r, f1, f3);
    EXPECT_TRUE(ors1.doesOverlap(f2, f2));
    EXPECT_FALSE(ors1.doesOverlap(f4, f4));
    RingSection ors2 = buildORS(r, f1, f4);
    EXPECT_TRUE(ors2.doesOverlap(f2, f3));
    RingSection ors3 = buildORS(r, f2, f3);
    EXPECT_TRUE(ors3.doesOverlap(f2, f3));
    EXPECT_FALSE(ors3.doesOverlap(f1, f1));
    EXPECT_FALSE(ors3.doesOverlap(f4, f4));
}

TEST(RingSectionTest, TestGlobalOffset1) {
    // Try a global offset having same direction and position as the 
    // start position and direction. following handleOffset, the offset should
    // be a nullOp
    Offset testOffset = Offset::globalCartesianOffset("aname", 
                                                      Vector_t(1.0, 2.0, 3.0),
                                                      Vector_t(4.0, 5.0, 6.0));
    EXPECT_FALSE(testOffset.getIsLocal());
    RingSection section;
    //handleOffset needs start position/normal to be set
    section.setStartPosition(Vector_t(1.0, 2.0, 3.0));
    section.setStartNormal(Vector_t(4.0, 5.0, 6.0));
    section.setComponent(&testOffset); // this borrows the offset pointer
    try {
        section.handleOffset();
    } catch (ClassicException& exc) {
        EXPECT_TRUE(false) << "threw an exception: " << exc.what();
    }
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(testOffset.getEndPosition()[i], Vector_t(0, 0, 0)[i]) << i;
        EXPECT_EQ(testOffset.getEndDirection()[i], Vector_t(0, 1, 0)[i]) << i;
    }
    EXPECT_TRUE(testOffset.getIsLocal());
}

TEST(RingSectionTest, TestGlobalOffset2) {
    Vector_t startPos(1.0, 2.0, 3.0);
    Vector_t endPos(2.0, 3.0, 3.0);
    Vector_t startDir(-1.0, 0.0, 0.0);
    Vector_t endDir(-1.0, -1.0, 0.0);
    Offset testOffset =
        Offset::globalCartesianOffset("aname", endPos, endDir);
    RingSection section;
    //handleOffset needs start position/normal to be set
    section.setStartPosition(startPos);
    section.setStartNormal(startDir);
    section.setComponent(&testOffset); // this borrows the offset pointer
    section.handleOffset();
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(testOffset.getEndPosition()[i], Vector_t(1, -1, 0)[i], 1e-6) << i;
        EXPECT_NEAR(testOffset.getEndDirection()[i],
                  Vector_t(-1/std::sqrt(2), 1/std::sqrt(2), 0)[i], 1e-6) << i;
    }
    EXPECT_TRUE(testOffset.getIsLocal());
}
