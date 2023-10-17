//
// Unit tests for class Ring
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
#include "gtest/gtest.h"

#include "Algorithms/PartData.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Ring.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

#include "opal_src/Utilities/MockComponent.h"
#include "opal_test_utilities/SilenceTest.h"

#include <iostream>
#include <sstream>

// generate a set of weird, but closed, elements
// reaches theta sum after 16 elements
class OffsetFactory {
  public:
    OffsetFactory(double radius=1., int start=0, double thetaSum=-1.) {
        i_m = start;
        radius_m = radius;
        thetaSum_m = thetaSum;
        if (thetaSum_m < 0.)
            thetaSum_m = 2.*Physics::pi;
        nextIsMock_m = true;
    }

    // generate a set of weird, but closed, Offset elements
    // reaches theta sum after 16 elements
    Component* yield() {
        int cell = i_m % 8+1;
        double theta = thetaSum_m*cell/36./2.;
        double length = 2.*sin(theta/2.)*radius_m;
        i_m++;
        Offset* off = new Offset(Offset::localCylindricalOffset("offset1", 0., theta, length));
        offVec_m.push_back(off);
        return off;
    }

    // generate a set of weird, but closed, elements
    // alternate MockComponent (generates a field, straight 2 mm long) with
    // offset components from yield()
    // reaches theta sum after 12 elements
    Component* yieldComp1() {
        nextIsMock_m = !nextIsMock_m;
        if (nextIsMock_m) {
            Offset* off = new Offset(Offset::localCylindricalOffset("offset2", Physics::pi/6., Physics::pi/6., 2.));
            offVec_m.push_back(off);
            return off;
        }
        Offset* off = new Offset(Offset::localCylindricalOffset("offset3", 0., 0., 1));
        offVec_m.push_back(off);
        MockComponent* mock = new MockComponent();
        mock->geom_m = &off->getGeometry();
        mockVec_m.push_back(mock);
        return mock;
    }

    // generate a set of sector magnets; the geometry is defined so that they are
    // exact sector magnets
    Component* yieldComp2() {
        double f = Physics::pi/20.;
        Offset* off = new Offset(
           Offset::localCylindricalOffset("offset4", f, f, 2.*radius_m*sin(f)));
        offVec_m.push_back(off);
        MockComponent* mock = new MockComponent();
        mock->geom_m = &off->getGeometry();
        mockVec_m.push_back(mock);
        return mock;
    }

    int i_m;
    bool nextIsMock_m;
    double radius_m;
    double thetaSum_m;
    // keep the offset alive for the life of the test
    std::vector<Offset*> offVec_m;
    std::vector<MockComponent*> mockVec_m;
};

TEST(RingTest, TestConstructDestruct) {
    // something here? someday...
}

TEST(RingTest, TestAppend1) {
    OpalTestUtilities::SilenceTest silencer;

    try {
        double radius = 5.;
        Ring ring("my_ring");
        ring.setLatticeRInit(radius);
        ring.setLatticePhiInit(Physics::pi/2.);
        ring.setLatticeThetaInit(0.);
        ring.setSymmetry(1);
        ring.setIsClosed(true);
        Offset off = Offset::localCylindricalOffset("cyl1", 0., Physics::pi/6., 1);
        ring.appendElement(off);
        // note phiInit is at pi/2 i.e. on the y axis pointing towards negative x
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i), 
                        Vector_t(-1, 5., 0.)(i),
                        1e-6) << i;
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(-cos(Physics::pi/6.),
                                                          -sin(Physics::pi/6.),
                                                          0.)(i), 1e-6) << i;
        }
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i),
                        Vector_t(-1-cos(Physics::pi/6.),
                                  5-sin(Physics::pi/6.), 0.)(i), 1e-6) << i;
            Vector_t expected(-cos(Physics::pi/3.), -sin(Physics::pi/3.), 0.);
            EXPECT_NEAR(ring.getNextNormal()(i),
                        expected(i),
                        1e-6) << i;
        }
    } catch (OpalException& exc) {
        std::cerr << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n";
    }
}

TEST(RingTest, TestAppend2) {
    OpalTestUtilities::SilenceTest silencer;

    try {
        double radius = 5.;
        Ring ring("my_ring");
        ring.setLatticeRInit(radius);
        ring.setLatticePhiInit(0.);
        ring.setLatticeThetaInit(0.);
        ring.setSymmetry(1);
        ring.setIsClosed(true);
        Offset off = Offset::localCylindricalOffset("cyl1", Physics::pi/24., Physics::pi/8., 1);
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i),
                        Vector_t(5.-sin(Physics::pi/24.),
                                 cos(Physics::pi/24.), 0.)(i), 1e-6)
                << i << "\n";
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(-sin(Physics::pi/6.),
                                                           cos(Physics::pi/6.),
                                                           0.)(i), 1e-6)
                << i << "\n";
        }
        ring.appendElement(off);
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(-1., 0., 0.)(i), 1e-6)
                << i << "\n";
        }
    } catch (OpalException& exc) {
        std::cerr << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n";
    }
}

TEST(RingTest, TestAppend3) {
    OpalTestUtilities::SilenceTest silencer;

    try {
        double radius = 5.;
        Ring ring("my_ring");
        ring.setLatticeRInit(radius);
        ring.setLatticePhiInit(0.);
        ring.setLatticeThetaInit(0.);
        ring.setSymmetry(1);
        ring.setIsClosed(true);
        Offset off = Offset::localCylindricalOffset("cyl1", 0., Physics::pi/6., 1.);
        for (size_t i = 0; i < 12; ++i) {
            ring.appendElement(off);
        }
        ring.lockRing();
    } catch (OpalException& exc) {
        std::cerr << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n";
    }
}

TEST(RingTest, TestLatticeRInitPhiInit) {
    OpalTestUtilities::SilenceTest silencer;

    for (double phi = -2.*Physics::pi;
         phi < 2.*Physics::pi;
         phi += Physics::pi/6.) {
        for (double theta = -2.*Physics::pi;
             theta < 2.*Physics::pi;
             theta += Physics::pi/6.) {
            for (double radius = 1.; radius < 5.; radius += 1.) {
                Ring ring("my_ring");
                ring.setLatticeRInit(radius);
                ring.setLatticePhiInit(phi);
                ring.setLatticeThetaInit(theta);
                Vector_t pos = ring.getNextPosition();
                Vector_t refPos(radius*cos(phi), radius*sin(phi), 0.);
                for (size_t i = 0; i < 3; ++i) {
                    EXPECT_EQ(pos(i), refPos(i))
                        << i << " f: " << phi
                        << " t: " << theta << " r: " << radius << "\n";
                }
                Vector_t norm = ring.getNextNormal();
                Vector_t refNorm(-sin(phi+theta), cos(phi+theta), 0.);
                for (size_t i = 0; i < 3; ++i) {
                    EXPECT_EQ(norm(i), refNorm(i))
                        << i << " f: " << phi
                        << " t: " << theta << " r: " << radius << "\n";
                }
            }
        }
    }
}


// Check that we get the field lookup correct accounting for the position of the
// field element
TEST(RingTest, TestApply) {
    OpalTestUtilities::SilenceTest silencer;

    for (double phi = 0.0; phi < 3*Physics::pi+1e-3; phi += Physics::pi/3.0) {
        Ring ring("my_ring");
        try {
            double radius = 2.*Units::mm2m;
            PartData data;
            PartBunch bunch(&data);
            ring.setRefPartBunch(&bunch);
            ring.setLatticeRInit(radius);
            ring.setLatticePhiInit(phi);
            ring.setLatticeThetaInit(0.);
            ring.setSymmetry(1);
            ring.setIsClosed(false);
            OffsetFactory fac(radius);
            ring.appendElement(*fac.yieldComp1());
            ring.lockRing();

            // position should always be 0.5, -0.75, 0.25 in the local
            // (OPAL-T) coordinate system. B Field should be the vector between
            // the field looku position and the start position of the element
            Vector_t position(2.5*Units::mm2m, 0.25*Units::mm2m, -0.75*Units::mm2m);
            Vector_t rotated_position, P, B, E;
            rotated_position[0] = position[0]*cos(phi) - position[1]*sin(phi);
            rotated_position[1] = position[0]*sin(phi) + position[1]*cos(phi);
            rotated_position[2] = position[2];
            bool outOfBounds = ring.apply(rotated_position, P, 0.0,  E, B);
            Vector_t refB = rotated_position - (ring.getSection(0)->getStartPosition());
            for (size_t i = 0; i < 3; ++i) {
                EXPECT_FALSE(outOfBounds);
                EXPECT_NEAR(B[i], refB[i], 1e-6) << phi << " " << i;
                EXPECT_NEAR(E[i], -refB[i], 1e-6) << phi << " " << i;
            }
        } catch (ClassicException& exc) {
            std::cout << exc.what() << std::endl;
            EXPECT_TRUE(false) << "Threw an exception\n";
        }
    }
}
