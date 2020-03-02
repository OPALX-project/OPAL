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

#include "Algorithms/PartData.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/Offset.h"
#include "opal_src/Utilities/MockComponent.h"
#include "AbsBeamline/Ring.h"
#include "Utilities/OpalException.h"

#include "opal_test_utilities/SilenceTest.h"

#include <iostream>
#include <sstream>

const double units=1e3;

// generate a set of weird, but closed, elements
// reaches theta sum after 16 elements
class OffsetFactory {
  public:
    OffsetFactory(double radius=1., int start=0, double thetaSum=-1.) {
        i_m = start;
        radius_m = radius/units;
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
            Offset* off = new Offset(Offset::localCylindricalOffset("offset2", Physics::pi/6., Physics::pi/6., 2./units));
            offVec_m.push_back(off);
            return off;
        }
        Offset* off = new Offset(Offset::localCylindricalOffset("offset3", 0., 0., 1./units));
        offVec_m.push_back(off);
        MockComponent* mock = new MockComponent();
        mock->geom_m = &off->getGeometry();
        mockVec_m.push_back(mock);
        return mock;
    }

    // generate a set of sector magnets; the geometry is defined so that they are
    // exact sector magnets
    Component* yieldComp2() {
        std::cerr << "YIELDCOMP2" << std::endl;
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
        Offset off = Offset::localCylindricalOffset("cyl1", 0., Physics::pi/6., 1./units);
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i), Vector_t(5., -1., 0.)(i), 1e-6);
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(-sin(Physics::pi/6.),
                                                          -cos(Physics::pi/6.),
                                                          0.)(i), 1e-6);
        }
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i),
                        Vector_t(5.-sin(Physics::pi/6.),
                                 -1.-cos(Physics::pi/6.), 0.)(i), 1e-6);
            Vector_t expected(-sin(Physics::pi/3.), -cos(Physics::pi/3.), 0.);
            EXPECT_NEAR(ring.getNextNormal()(i),
                        expected(i),
                        1e-6);
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
        Offset off = Offset::localCylindricalOffset("cyl1", Physics::pi/24., Physics::pi/8., 1./units);
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextPosition()(i),
                        Vector_t(cos(Physics::pi/24.),
                                 5.-sin(Physics::pi/24.), 0.)(i), 1e-6)
                << i << "\n";
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(cos(Physics::pi/6.),
                                                         -sin(Physics::pi/6.),
                                                          0.)(i), 1e-6)
                << i << "\n";
        }
        ring.appendElement(off);
        ring.appendElement(off);
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(ring.getNextNormal()(i), Vector_t(0., -1., 0.)(i), 1e-6)
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
                Vector_t refPos(radius*sin(phi), radius*cos(phi), 0.);
                for (size_t i = 0; i < 3; ++i) {
                    EXPECT_EQ(pos(i), refPos(i))
                        << i << " f: " << phi
                        << " t: " << theta << " r: " << radius << "\n";
                }
                Vector_t norm = ring.getNextNormal();
                Vector_t refNorm(cos(phi+theta), -sin(phi+theta), 0.);
                for (size_t i = 0; i < 3; ++i) {
                    EXPECT_EQ(norm(i), refNorm(i))
                        << i << " f: " << phi
                        << " t: " << theta << " r: " << radius << "\n";
                }
            }
        }
    }
}

// Check that we get the bounding box and rotation correct
TEST(RingTest, TestApply) {
    OpalTestUtilities::SilenceTest silencer;

    Ring ring("my_ring");
    double metres = 1e-3;
    try {
        double radius = 2.*(2.*sin(Physics::pi/6.)+1.*sin(Physics::pi/3.)+1.0);
        PartData data;
        PartBunch bunch(&data);
        ring.setRefPartBunch(&bunch);
        ring.setLatticeRInit(radius-2.);
        ring.setLatticePhiInit(0.);
        ring.setLatticeThetaInit(0.);
        ring.setSymmetry(1);
        ring.setIsClosed(true);
        OffsetFactory fac(radius);
        for (size_t i = 0; i < 12; ++i) {
            ring.appendElement(*fac.yieldComp1());
        }
        ring.lockRing();
        // check that we get a MockComponent rotated thru 180 degrees
        for (double x = -1.0001; x < 2.; x += 0.1) {
            double y = -2.2;
            Vector_t pos(x, y, -0.5);
            pos *= metres;
            Vector_t zero(0.0);
            Vector_t centroid(0., 0., 0), B(0., 0., 0), E(0., 0., 0);
            double t = 0;
            std::cout << "Apply pos:" << pos << std::endl;
            EXPECT_FALSE(ring.apply(pos, zero, t, E, B));
            std::cout << "Yields B: " << B << " E: " << E << std::endl;
            Vector_t BRef(0.0, 0.0, 0.0);
            if (x > 0. and x < 1.)
                BRef = Vector_t(x-1., y+2., -0.5);
            std::cout << "Expected B: " << BRef << std::endl;
            for (int i = 0; i < 3; ++i) {
                EXPECT_NEAR(B(i), BRef(i), 1e-6)
                               << "component " << i << " for pos " << pos;
                EXPECT_NEAR(E(i), -BRef(i), 1e-6)
                               << "component " << i << " for pos " << pos;
            }
        }
        // check that we get something reasonable for all phi
        for (double phi = 0.; phi < 2.*Physics::pi+0.1; phi += Physics::pi/100.) {
            Vector_t pos(radius/2.*sin(phi), radius/2.+radius/2.*cos(phi), 0.5);
            Vector_t centroid, B, E;
            EXPECT_FALSE(ring.apply(pos, Vector_t(0.0), 0., E, B)); // check we don't throw for all angles
            // std::cout << phi << " " << pos << " " << B << std::endl;
        }
    } catch (OpalException& exc) {
        std::cout << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n";
    }
}

// Check that we get the bounding box correct - for exact sector geometry
TEST(RingTest, TestApply2) {
    OpalTestUtilities::SilenceTest silencer;

    Ring ring("my_ring");
    double metres = 1e-3;
    try {
        double radius = 1.5;
        PartData data;
        PartBunch bunch(&data);
        ring.setRefPartBunch(&bunch);
        ring.setLatticeRInit(radius);
        ring.setLatticePhiInit(7.*Physics::pi/4.);
        ring.setLatticeThetaInit(0.);
        ring.setSymmetry(1);
        ring.setIsClosed(true);
        OffsetFactory fac(radius);
        for (size_t i = 0; i < 20; ++i) {
            ring.appendElement(*fac.yieldComp2());
        }
        ring.lockRing();
        for (double phi = 0.001; phi < 2.*Physics::pi+0.1; phi += Physics::pi/50.) {
            Vector_t pos((radius+0.5)*sin(phi), (radius+0.5)*cos(phi), -0.5);
            pos *= metres;
            Vector_t centroid, B, E;
            std::vector<RingSection*> sections = ring.getSectionsAt(pos);
            EXPECT_FALSE(ring.apply(pos, Vector_t(0.0), 0., E, B));
            // check we don't throw for all angles
            // a few are coming out with Bz = 1. instead of Bz = 0.5; looks like
            // floating point precision issue? It's okay, Ring is not
            // responsible for bounding the field, Components are.
            EXPECT_GE(-B(2), 0.1);
            EXPECT_LE(-B(2), 1.1);
            // std::cout << phi << " " << pos << " " << B << " " << sections.size() << std::endl;
        }
    } catch (OpalException& exc) {
        std::cout << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n" ;
    }
    // Now apply symmetry 2x10 fields instead of 20x1
    Ring ring2("my_ring");
    try {
        double radius = 1.5;
        PartData data;
        PartBunch bunch(&data);
        ring2.setRefPartBunch(&bunch);
        ring2.setLatticeRInit(radius);
        ring2.setLatticePhiInit(7.*Physics::pi/4.);
        ring2.setLatticeThetaInit(0.);
        ring2.setSymmetry(10);
        ring2.setIsClosed(true);
        OffsetFactory fac(radius);
        for (size_t i = 0; i < 2; ++i) {
            ring2.appendElement(*fac.yieldComp2());
        }
        ring2.lockRing();
        for (double phi = 0.001; phi < 2.*Physics::pi+0.1; phi += Physics::pi/50.) {
            Vector_t pos((radius+0.5)*sin(phi), (radius+0.5)*cos(phi), 0.5);
            Vector_t centroid, B1, B2, E;
            std::vector<RingSection*> sections = ring2.getSectionsAt(pos);
            ring.apply(pos, Vector_t(0.0), 0., E, B1);
            ring2.apply(pos, Vector_t(0.0), 0., E, B2);
            EXPECT_NEAR(B1(2), B2(2), 1e-6);
            // std::cout << phi << " " << pos << " " << B << " " << sections.size() << std::endl;
        }
    } catch (OpalException& exc) {
        std::cout << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n" ;
    }
    // Now overlapping - we have two elements in each position, should get twice
    // the field
    Ring ring3("my_ring");
    try {
        double radius = 1.5;
        PartData data;
        PartBunch bunch(&data);
        ring3.setRefPartBunch(&bunch);
        ring3.setLatticeRInit(radius);
        ring3.setLatticePhiInit(7.*Physics::pi/4.);
        ring3.setLatticeThetaInit(0.);
        ring3.setIsClosed(true);
        OffsetFactory fac(radius);
        for (size_t i = 0; i < 40; ++i) {
            ring3.appendElement(*fac.yieldComp2());
        }
        ring3.lockRing();
        for (double phi = 0.001; phi < 2.*Physics::pi+0.1; phi += Physics::pi/50.) {
            Vector_t pos((radius+0.5)*sin(phi), (radius+0.5)*cos(phi), 0.5);
            Vector_t centroid, B1, B2, E;
            std::vector<RingSection*> sections = ring3.getSectionsAt(pos);
            ring.apply(pos, Vector_t(0.0), 0., E, B1);
            ring3.apply(pos, Vector_t(0.0), 0., E, B2);
            EXPECT_NEAR(2.*B1(2), B2(2), 1e-6);
            // std::cout << phi << " " << pos << " " << B << " " << sections.size() << std::endl;
        }
    } catch (OpalException& exc) {
        std::cout << exc.what() << std::endl;
        EXPECT_TRUE(false) << "Threw an exception\n";
    }
}

void testField(double s, double r, double y, double phi,
               double bx, double by, double bz, double /*tol*/) {
    double radius = 2.;
    Ring ring("test");
    PartData data;
    PartBunch bunch(&data);
    ring.setRefPartBunch(&bunch);
    ring.setLatticeRInit(radius);
    ring.setLatticePhiInit(phi);
    ring.setLatticeThetaInit(0.);
    ring.setSymmetry(1);
    ring.setIsClosed(false);
    MockComponent2 mock;
    Offset off = Offset(Offset::localCylindricalOffset("offset", 0., 0., 10.));
    mock.geom_m = &off.getGeometry();
    ring.appendElement(mock);
    ring.lockRing();
    Vector_t centroid, E, B;
    Vector_t pos(radius*sin(phi)+s*cos(phi)+r*sin(phi),
                 radius*cos(phi)-s*sin(phi)+r*cos(phi),
                 y);
    pos *= 1e-3; // metres
    ring.apply(pos, Vector_t(0.0), 0., E, B);
    EXPECT_NEAR(B(0), bx, 1e-6);
    EXPECT_NEAR(B(1), by, 1e-6);
    EXPECT_NEAR(B(2), bz, 1e-6);
    // std::cout << pos << " ** " << B << " ** " << Vector_t(bx, by, bz) << std::endl;
}

TEST(RingTest, TestApply3) {
    OpalTestUtilities::SilenceTest silencer;

    testField(0.1, 0., 0.2, 0., 3., 1., 2., 1e-6);
    testField(0.1, 0., 0.2, Physics::pi, -3., -1., 2., 1e-6);
    testField(0.1, 0., 0.2, Physics::pi/2., 1., -3., 2., 1e-6);
    testField(0.1, 0., 0.2, 3.*Physics::pi/2., -1., 3., 2., 1e-6);
    testField(0.1, 0.15, 0.2, Physics::pi/6.,
              3.*cos(Physics::pi/6)+1.*sin(Physics::pi/6),
              -3.*sin(Physics::pi/6)+1.*cos(Physics::pi/6), 2., 1e-6);
}