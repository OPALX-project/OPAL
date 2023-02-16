#include "gtest/gtest.h"

#include "TrimCoils/TrimCoilBFit.h"
#include "TrimCoils/TrimCoilPhaseFit.h"
#include "TrimCoils/TrimCoilMirrored.h"

#include "Physics/Units.h"
#include "Utilities/Util.h"

#include "opal_test_utilities/SilenceTest.h"

const double margin = 1e-7;

// Test TrimCoilBFit on zeros
TEST(TrimCoil, TrimCoilBFitZeros)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 0.0;
    double rmin = 1000; // mm
    double rmax = 2000;
    double phi  = 0.0;  // rad

    TrimCoilBFit myTrimCoil(bmax, rmin, rmax, {}, {}, {}, {});

    const double one = 1.0;
    double br = one, bz = one;
    myTrimCoil.applyField((rmin + rmax) * Units::mm2m / 2.0, 1, 0, &br, &bz);
    // not changed since bmax 0.0
    EXPECT_NEAR(br, one, margin);
    EXPECT_NEAR(bz, one, margin);

    bmax = 1.0;
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {}, {}, {}, {});

    myTrimCoil.applyField(rmin * Units::mm2m - 1, 1, phi, &br, &bz);
    // not changed since r outside range
    EXPECT_NEAR(br, one, margin);
    EXPECT_NEAR(bz, one, margin);

    myTrimCoil.applyField(rmax * Units::mm2m + 1, 1, phi, &br, &bz);
    // not changed since r outside range
    EXPECT_NEAR(br, one, margin);
    EXPECT_NEAR(bz, one, margin);

    myTrimCoil.applyField(rmax * Units::mm2m - 1, 1, phi, &br, &bz);
    // default constant field
    EXPECT_NEAR(br, one, margin);
    EXPECT_NEAR(bz, 11.0, margin); // 1 + 10*1
}

// TrimCoilBFit class
TEST(TrimCoil, TrimCoilBFit)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;
    double phi  = 0.0;  // rad

    // polynom 1 + 2*x + 3*x^2
    TrimCoilBFit myTrimCoil(bmax, rmin, rmax, {1.0, 2.0, 3.0}, {}, {}, {});

    const double brStart = 1.0, bzStart = 1.0, zStart = 2.0;
    double br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, phi, &br, &bz);
    // 1 + 10*(1 + 2*2 + 3*2*2) = 171 ; 1 + 10*2*(2 + 6*2) = 281
    EXPECT_NEAR(bz, 171.0, margin);
    EXPECT_NEAR(br, 281.0, margin);

    // rational function (4 + 3*x) / (1 + 2*x)
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0}, {}, {});

    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, phi, &br, &bz);

    double bzSolution = 21.0; // = 1 + 10*(4+3*2) / (1+2*2)
    double brSolution = -3.0; // = 1 + 10*2*(-0.2)

    EXPECT_NEAR(bz, bzSolution, margin);
    EXPECT_NEAR(br, brSolution, margin);

    phi = Util::angle_0to2pi(1.0);
    EXPECT_NEAR(phi, 1, margin);
    phi = Util::angle_0to2pi(-1.0);
    EXPECT_NEAR(phi, -1+Physics::two_pi, margin);

    // test phi angles
    myTrimCoil.setAzimuth(10,180);
    bool valid = Util::angleBetweenAngles(0.0, 10*Units::deg2rad, 180*Units::deg2rad);
    EXPECT_FALSE(valid);
    valid = Util::angleBetweenAngles(1.0, 10*Units::deg2rad, 180*Units::deg2rad);
    EXPECT_TRUE(valid);
    valid = Util::angleBetweenAngles(3.2, 10*Units::deg2rad, 180*Units::deg2rad);
    EXPECT_FALSE(valid);
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, 3.2, &br, &bz); // outside range
    EXPECT_NEAR(bz, bzStart, margin);
    myTrimCoil.applyField(2.0, zStart, 0.0, &br, &bz); // outside range
    EXPECT_NEAR(bz, bzStart, margin);
    myTrimCoil.applyField(2.0, zStart, 1.0, &br, &bz); //  inside range
    EXPECT_NEAR(bz, bzSolution, margin);

    // test phi angles: first larger
    myTrimCoil.setAzimuth(180,20);
    valid = Util::angleBetweenAngles(0.0, 180*Units::deg2rad, 20*Units::deg2rad);
    EXPECT_TRUE(valid);
    valid = Util::angleBetweenAngles(1.0, 180*Units::deg2rad, 20*Units::deg2rad);
    EXPECT_FALSE(valid);
    valid = Util::angleBetweenAngles(3.2, 180*Units::deg2rad, 20*Units::deg2rad);
    EXPECT_TRUE(valid);
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, 1.0, &br, &bz); // outside range
    EXPECT_NEAR(bz, bzStart, margin);
    myTrimCoil.applyField(2.0, zStart, 3.2, &br, &bz); //  inside range
    EXPECT_NEAR(bz, bzSolution, margin);
    myTrimCoil.applyField(2.0, zStart, 0.0, &br, &bz); //  inside range
    EXPECT_NEAR(bz, 2*bzSolution - bzStart, margin);

    // test phi angles: negative
    myTrimCoil.setAzimuth(10,-180);
    double pi = Util::angle_0to2pi(-180*Units::deg2rad);
    valid = Util::angleBetweenAngles(0.0, 10*Units::deg2rad, pi);
    EXPECT_FALSE(valid);
    valid = Util::angleBetweenAngles(1.0, 10*Units::deg2rad, pi);
    EXPECT_TRUE(valid);
    valid = Util::angleBetweenAngles(3.2, 10*Units::deg2rad, pi);
    EXPECT_FALSE(valid);
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, 3.2, &br, &bz); // outside range
    EXPECT_NEAR(bz, bzStart, margin);
    myTrimCoil.applyField(2.0, zStart, 0.0, &br, &bz); // outside range
    EXPECT_NEAR(bz, bzStart, margin);
    myTrimCoil.applyField(2.0, zStart, 1.0, &br, &bz); //  inside range
    EXPECT_NEAR(bz, bzSolution, margin);

    // Test phi function
    // same rational function
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {}, {}, {4.0, 3.0}, {1.0, 2.0});
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(1.0, zStart, 2.0, &br, &bz);
    EXPECT_NEAR(bz, bzSolution, margin);
    EXPECT_NEAR(br, brSolution, margin);
    // Both r and phi
    myTrimCoil = TrimCoilBFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0}, {4.0, 3.0}, {1.0, 2.0});
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, 2.0, &br, &bz);
    EXPECT_NEAR(bz,  41.0, margin); // 1 + 10 * ((4+3*2) / (1+2*2)) **2
    EXPECT_NEAR(br, -15.0, margin);
}

// TrimCoilPhaseFit class
TEST(TrimCoil, TrimCoilPhaseFit)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;
    double phi  = 0.0;

    // polynom 1 + 2*x + 3*x^2
    TrimCoilPhaseFit myTrimCoil(bmax, rmin, rmax, {1.0, 2.0, 3.0}, {}, {}, {});
    const double brStart = 1.0, bzStart = 1.0, zStart = 2.0;
    double br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, phi, &br, &bz);
    // 1 - 10*(2 + 6*2) = -139; 1 - 2*10*6 = -119
    EXPECT_NEAR(bz, -139.0, margin);
    EXPECT_NEAR(br, -119.0, margin);

    // rational function (4 + 3*x) / (1 + 2*x)
    myTrimCoil = TrimCoilPhaseFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0}, {}, {});

    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, phi, &br, &bz);
    double bzSolution =  3.0; // = 1 - 10*(-0.2)
    double brSolution = -2.2; // = 1 + 2*10*2*(-0.2*2) / 5
    EXPECT_NEAR(bz, bzSolution, margin);
    EXPECT_NEAR(br, brSolution, margin);
    // Test Phi
    // same rational function
    myTrimCoil = TrimCoilPhaseFit(bmax, rmin, rmax, {}, {}, {4.0, 3.0}, {1.0, 2.0});
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(1.0, zStart, 2.0, &br, &bz);
    EXPECT_NEAR(bz, bzSolution, margin);
    EXPECT_NEAR(br, brSolution, margin);
    // Both r and phi
    myTrimCoil = TrimCoilPhaseFit(bmax, rmin, rmax, {4.0, 3.0}, {1.0, 2.0}, {4.0, 3.0}, {1.0, 2.0});
    br = brStart, bz = bzStart;
    myTrimCoil.applyField(2.0, zStart, 2.0, &br, &bz);
    EXPECT_NEAR(bz,   9.0, margin);
    EXPECT_NEAR(br, -13.4, margin);
}

// TrimCoilMirrored class
TEST(TrimCoil, TrimCoilMirrored)
{
    OpalTestUtilities::SilenceTest silencer;

    double bmax = 1.0; // 1 T will be converted to 10 kG
    double rmin = 0.0;
    double rmax = 3000.0;
    double phi  = 0.0;
    double bslope = 1./6.;

    TrimCoilMirrored myTrimCoil(bmax, rmin, rmax, bslope);

    double br = 1.0, bz = 1.0;
    myTrimCoil.applyField(2.0, 2.0, phi, &br, &bz);

    EXPECT_NEAR(bz,-6.1943868603626751, 1e-6);
    EXPECT_NEAR(br, 1.0032755233321968, 1e-6);
}

