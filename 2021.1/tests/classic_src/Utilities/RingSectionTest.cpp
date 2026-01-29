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

#include <algorithm>

#include "gtest/gtest.h"
#include "opal_src/Utilities/MockComponent.h"
#include "Physics/Physics.h"
#include "Utilities/RingSection.h"

#include "opal_test_utilities/SilenceTest.h"

TEST(RingSectionTest, TestConstructDestruct) {
    OpalTestUtilities::SilenceTest silencer;

    RingSection ors;
    MockComponent* compNull = NULL;
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
    Vector_t centre(-1.33, +1.66, 0.);
    for (double theta = -3.*Physics::pi; theta < 3.*Physics::pi; theta += Physics::pi/6.) {
        Vector_t orientation(0., 0., theta);
        ors.setComponentOrientation(orientation);
        ors.setComponentPosition(centre);
        double c = cos(orientation(2));
        double s = -sin(orientation(2));
        for (double x = 0.01; x < 1.; x += 0.1)
            for (double y = 0.01; y < 1.; y += 0.1)
                for (double z = -0.01; z > -1.; z -= 0.1) {
                    Vector_t offset(c*x+s*y, -s*x+c*y, z);
                    Vector_t pos = centre+offset;
                    Vector_t centroid, B, E;
                    double t = 0;
                    EXPECT_FALSE(ors.getFieldValue(pos, centroid, t, E, B));
                    Vector_t bfield(c*x+s*y, -s*x+c*y, z);
                    for (int l = 0; l < 3; ++l) {
                        EXPECT_NEAR(B(l), +bfield(l), 1e-6);
                        EXPECT_NEAR(E(l), -bfield(l), 1e-6);
                    }
                }
    }
}

namespace {
    bool sort_comparator(Vector_t v1, Vector_t v2) {
        if (fabs(v1(0) - v2(0)) < 1e-6) {
            if (fabs(v1(1) - v2(1)) < 1e-6) {
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
    bbRef.push_back(Vector_t(0.99*sqrt(10)/(-sqrt(17))+3.,
                             0.99*sqrt(10)*4./(+sqrt(17))-1., 99.));
    bbRef.push_back(Vector_t(0.99*sqrt(10)/(+sqrt(17))+3.,
                             0.99*sqrt(10)*4./(-sqrt(17))-1., 99.));
    bbRef.push_back(Vector_t(0.99*sqrt(5)/(+sqrt(2))+2.,
                             0.99*sqrt(5)/(+sqrt(2))+1., 77.));
    bbRef.push_back(Vector_t(0.99*sqrt(5)/(-sqrt(2))+2.,
                             0.99*sqrt(5)/(-sqrt(2))+1., 77.));
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
        ors.setStartPosition(Vector_t(sin(phi1)*r, cos(phi1)*r, 0.));
        ors.setStartNormal(Vector_t(cos(phi1), -sin(phi1), 0.));
        ors.setEndPosition(Vector_t(sin(phi2)*r, cos(phi2)*r, 0.));
        ors.setEndNormal(Vector_t(cos(phi2), -sin(phi2), 0.));
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
