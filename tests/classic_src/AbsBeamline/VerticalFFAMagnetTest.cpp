//
// Unit tests for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//


#include "AbsBeamline/VerticalFFAMagnet.h"

#include <cmath>
#include <fstream>
#include <sstream>

#include "gtest/gtest.h"
#include "opal_test_utilities/SilenceTest.h"
#include "opal_test_utilities/Maxwell.h"


#include "Classic/AbsBeamline/EndFieldModel/Tanh.h"
#include "Classic/AbsBeamline/VerticalFFAMagnet.h"

class VerticalFFAMagnetTest : public ::testing::Test {
public:
    VerticalFFAMagnetTest() {
    }

    void SetUp( ) {
        magnet_m.reset(new VerticalFFAMagnet("test"));
        tanh_m = new endfieldmodel::Tanh(length_m*mm, length_m*mm/4., 20);
        magnet_m->setEndField(tanh_m);
        magnet_m->setMaxOrder(maxOrder_m);
        magnet_m->setB0(b0_m);
        magnet_m->setFieldIndex(k_m);
        magnet_m->setNegativeVerticalExtent(length_m/4.);
        magnet_m->setPositiveVerticalExtent(length_m*2.);
        magnet_m->setBBLength(length_m*4.);
        magnet_m->setWidth(length_m);
        magnet_m->initialise();
    }

    void TearDown( ) {
    }

    ~VerticalFFAMagnetTest() {
    }

    std::unique_ptr<VerticalFFAMagnet> magnet_m;
    double length_m = 1.0;
    double k_m = 1.0; // m^-1
    size_t maxOrder_m = 2;
    double b0_m = 1.0;
    endfieldmodel::Tanh* tanh_m; // magnet_m owns this memory
    const double mm = 1e3;
    const double tesla = 10.;
    OpalTestUtilities::SilenceTest silencer_m;
};

TEST_F(VerticalFFAMagnetTest, ConstructorTest) {
    // check the constructor fills values okay
    EXPECT_EQ(magnet_m->getMaxOrder(), maxOrder_m);
    EXPECT_EQ(magnet_m->getEndField(), tanh_m); // borrowed pointer
    EXPECT_NEAR(magnet_m->getB0(), b0_m, 1e-6);
    EXPECT_NEAR(magnet_m->getFieldIndex(), k_m, 1e-6);
    EXPECT_NEAR(magnet_m->getNegativeVerticalExtent(), length_m/4, 1e-6);
    EXPECT_NEAR(magnet_m->getPositiveVerticalExtent(), length_m*2, 1e-6);
    EXPECT_NEAR(magnet_m->getBBLength(), length_m*4, 1e-6);
    EXPECT_NEAR(magnet_m->getWidth(), length_m, 1e-6);
}

TEST_F(VerticalFFAMagnetTest, MidplaneFieldTest) {
    // check the dipole field on the midplane meets vffa description
    // check the horizontal field is 0
    // The longitudinal field is generally non-zero, and is checked by maxwell
    // test.
    for (double z = -length_m*2.; z < length_m*2.; z += length_m/10.) {
        for (double y = -length_m/4.; y < length_m*2.; y += length_m/10.) {
            Vector_t position(0., y, z+length_m*2.); // z is relative to centre
            position *= mm;
            Vector_t bfield;
            bool outOfBounds = magnet_m->getFieldValue(position, bfield);
            double bRef = b0_m*tanh_m->function(z*mm, 0)*exp(k_m*y)*10.;
            EXPECT_FALSE(outOfBounds);
            EXPECT_NEAR(bfield(0), 0., 1e-12);
            EXPECT_NEAR(bfield(1), bRef, 1e-12);
        }
      }
}

TEST_F(VerticalFFAMagnetTest, BoundingBoxTest) {
    // check getFieldValue returns true if outside the bounding box
    std::vector< std::pair<Vector_t, bool> > tests;
    tests.insert(tests.begin(), {
        std::pair<Vector_t, bool>(Vector_t(0., 0., 1e-6), false),
        std::pair<Vector_t, bool>(Vector_t(0., 0., -1e-6), true),
        std::pair<Vector_t, bool>(Vector_t(0., 0., length_m*4.*mm-1e-6), false),
        std::pair<Vector_t, bool>(Vector_t(0., 0., length_m*4.*mm+1e-6), true),
        std::pair<Vector_t, bool>(Vector_t(0., -length_m/4.*mm+1e-6, length_m*2*mm), false),
        std::pair<Vector_t, bool>(Vector_t(0., -length_m/4.*mm-1e-6, length_m*2*mm), true),
        std::pair<Vector_t, bool>(Vector_t(0., length_m*2.*mm-1e-6, length_m*2*mm), false),
        std::pair<Vector_t, bool>(Vector_t(0., length_m*2.*mm+1e-6, length_m*2*mm), true),
        std::pair<Vector_t, bool>(Vector_t(-length_m/2*mm+1e-6, 0., length_m*2*mm), false),
        std::pair<Vector_t, bool>(Vector_t(-length_m/2*mm-1e-6, 0., length_m*2*mm), true),
        std::pair<Vector_t, bool>(Vector_t(+length_m/2*mm-1e-6, 0., length_m*2*mm), false),
        std::pair<Vector_t, bool>(Vector_t(+length_m/2*mm+1e-6, 0., length_m*2*mm), true),
    });

    Vector_t bfield;
    for (auto a_test = tests.begin(); a_test != tests.end(); ++a_test) {
        bool outOfBounds = magnet_m->getFieldValue(a_test->first, bfield);
        EXPECT_EQ(outOfBounds, a_test->second) << a_test->first;
    }
}

TEST_F(VerticalFFAMagnetTest, MaxwellTest) {
    // check the field is maxwellian. This is key test for the field model
    double x = length_m*0.1;
    double y = length_m*0.05;
    double z = length_m*0.8;
    VerticalFFAMagnet* magnet = dynamic_cast<VerticalFFAMagnet*>(magnet_m->clone());
    double dr = 1e-3;
    MaxwellTest maxTest(Vector_t(dr, dr, dr), 1., magnet);
    //maxTest.printHeading(std::cerr);
    double divOld = 1.e9;
    std::cerr << "Testing Maxwellianness" << std::endl;
    for (size_t i = 1; i < 12; i += 1) {
        magnet->setMaxOrder(i);
        magnet->initialise();
        Vector_t pos(x*mm, y*mm, z*mm);
        //maxTest.printLine(std::cerr, pos, 0.);
        double div = maxTest.divB(pos, 0.);
        double curl = euclidean_norm(maxTest.curlB(pos, 0.));
        EXPECT_LT(std::abs(div), std::abs(divOld)) << i;
        EXPECT_LT(curl, 1e-11);
        std::cerr << "Max Order: " << i << " |curlB|: " << curl
                  << " DivB: " << div << std::endl;
        if ((i/2)*2 != i) { // div only improves on even orders
            divOld = div;
        }
    }
}

TEST_F(VerticalFFAMagnetTest, CoefficientsTest) {
    // check the coefficient vector is filled; values must be correct if the
    // maxwell test passes so we don't test here.
    magnet_m->setMaxOrder(10);
    magnet_m->initialise();
    std::vector< std::vector<double> > coeffs = magnet_m->getDfCoefficients();
    EXPECT_EQ(coeffs.size(), 11u);
    for (size_t i = 0; i < coeffs.size(); i++) {
        std::vector<double> cv = coeffs[i];
        if ((i/2)*2 == i) { // even
            EXPECT_EQ(cv.size(), i+1u);
        } else { // odd
            EXPECT_EQ(cv.size(), 0u);
        }
    }
}