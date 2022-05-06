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

#include <memory>

#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"
#include "Physics/Physics.h"

#include "Algorithms/PolynomialTimeDependence.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"

class VariableRFCavityFringeFieldTest : public ::testing::Test {
public:
    VariableRFCavityFringeFieldTest() {
        cav1 = VariableRFCavityFringeField("bob");
        // centre, end, max_order
        shared = std::shared_ptr<endfieldmodel::EndFieldModel>
                                      (new endfieldmodel::Tanh(250., 50., 21));
        cav1.setCavityCentre(0.500);
        cav1.setEndField(shared);
        PolynomialTimeDependence* time =
                       new PolynomialTimeDependence(std::vector<double>(1, 1.));
        std::shared_ptr<AbstractTimeDependence> timePtr(time);
        cav1.setAmplitudeModel(timePtr);
        cav1.setPhaseModel(timePtr);
        cav1.setFrequencyModel(timePtr);
        cav1.setHeight(0.500);
        cav1.setWidth(0.200);
        cav1.setLength(1.000);
        cav1.initialiseCoefficients();

        cav2 = cav1;
        PolynomialTimeDependence* phase2 =
                       new PolynomialTimeDependence(std::vector<double>(1, 0.));
        PolynomialTimeDependence* freq2 =
                       new PolynomialTimeDependence(std::vector<double>(1, 1000.));
        std::shared_ptr<AbstractTimeDependence> phase2Ptr(phase2);
        std::shared_ptr<AbstractTimeDependence> freq2Ptr(freq2);
        cav2.setPhaseModel(phase2Ptr);
        cav2.setFrequencyModel(freq2Ptr);
    }

    VariableRFCavityFringeField cav1;
    VariableRFCavityFringeField cav2;
    std::shared_ptr<endfieldmodel::EndFieldModel> shared;
    OpalTestUtilities::SilenceTest silencer;
};

TEST_F(VariableRFCavityFringeFieldTest, TestConstructor) {
    std::cerr << "Test Ctor" << std::endl;
    VariableRFCavityFringeField cav("bob");
    EXPECT_FLOAT_EQ(cav.getCavityCentre(), 0.);

    EXPECT_FALSE(cav.getEndField());
    std::cerr << "Test Ctor 2" << std::endl;
}

TEST_F(VariableRFCavityFringeFieldTest, TestSetGet) {
    EXPECT_FLOAT_EQ(cav1.getCavityCentre(), 0.5); // mm
    EXPECT_FLOAT_EQ(cav1.getWidth(), 0.2); // metres
    EXPECT_EQ(&(*(cav1.getEndField())), &(*shared));
    long int count = shared.use_count();
    EXPECT_EQ(count, 3);
    PolynomialTimeDependence* time = nullptr;
    EXPECT_NE(&(*(cav1.getFrequencyModel())), time);
}

TEST_F(VariableRFCavityFringeFieldTest, TestClone) {
    VariableRFCavityFringeField* cav2 =
                       dynamic_cast<VariableRFCavityFringeField*>(cav1.clone());
    EXPECT_FLOAT_EQ(cav2->getCavityCentre(), 0.5);
    EXPECT_EQ(&(*(cav2->getEndField())), &(*shared));
    long int count = shared.use_count();
    EXPECT_EQ(count, 4);
}

TEST_F(VariableRFCavityFringeFieldTest, TestAssignment) {
    VariableRFCavityFringeField cav3 = cav1;
    EXPECT_FLOAT_EQ(cav3.getCavityCentre(), 0.5);
    EXPECT_EQ(&(*(cav3.getEndField())), &(*shared));
    long int count = shared.use_count();
    EXPECT_EQ(count, 4);
}

TEST_F(VariableRFCavityFringeFieldTest, TestCopyConstructor) {
    VariableRFCavityFringeField cav4(cav1);
    EXPECT_FLOAT_EQ(cav4.getCavityCentre(), 0.5);
    EXPECT_EQ(&(*(cav4.getEndField())), &(*shared));
    long int count = shared.use_count();
    EXPECT_EQ(count, 4);
}

TEST_F(VariableRFCavityFringeFieldTest, TestApplyBoundingBox) {
    Vector_t centroid(0., 0., 0.);
    Vector_t B(0., 0., 0.);
    Vector_t E(0., 0., 0.);
    std::vector<Vector_t> rVectorFalse = {
        Vector_t(-99.9, -249.9, 0.1),
        Vector_t(-99.9, 249.9, 999.9),
    };
    for (size_t i = 0; i < rVectorFalse.size(); ++i) {
        bool outOfBounds = cav1.apply(rVectorFalse[i], centroid, 0., B, E);
        EXPECT_FALSE(outOfBounds);
    }

      std::vector<Vector_t> rVectorTrue = {
        Vector_t(-100.1, 0., 500.),
        Vector_t(+100.1, 0., 500.),
        Vector_t(0., -250.1, 500.),
        Vector_t(0., +250.1, 500.),
        Vector_t(0., 0., -0.1),
        Vector_t(0., 0., 1000.1),
    };
    for (size_t i = 0; i < rVectorTrue.size(); ++i) {
        bool outOfBounds = cav1.apply(rVectorTrue[i], centroid, 0., B, E);
        EXPECT_TRUE(outOfBounds) << rVectorTrue[i];
    }
}

void testFieldLookup(VariableRFCavityFringeField& cav, Vector_t R, double t, Vector_t E, Vector_t B) {
    Vector_t centroid(0., 0., 0.);
    Vector_t Btest(0., 0., 0.);
    Vector_t Etest(0., 0., 0.);
    cav.apply(R, centroid, t, Etest, Btest);
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_FLOAT_EQ(E[i], Etest[i]) << "\nR:" << R << " t: " << t
        << "\nE expected: " << E << " " << " E meas: " << Etest
        << "\nB expected: " << B << " " << " B meas: " << Btest << std::endl;
        EXPECT_FLOAT_EQ(B[i], Btest[i]) << "\nR:" << R
        << "\nE expected: " << E << " " << " E meas: " << Etest
        << "\nB expected: " << B << " " << " B meas: " << Btest << std::endl;
    }
}

void partial(VariableRFCavityFringeField& cav, Vector_t pos, double t, double delta, int var, Vector_t& dE, Vector_t& dB) {
    bool verbose = false;
    Vector_t centroid(0., 0., 0.);
    Vector_t Bplus(0., 0., 0.);
    Vector_t Bminus(0., 0., 0.);
    Vector_t Eplus(0., 0., 0.);
    Vector_t Eminus(0., 0., 0.);
    Vector_t posPlus(pos), posMinus(pos);
    double tPlus(t), tMinus(t);
    if (var == 3) {
        tMinus -= delta;
        tPlus += delta;
    } else if (var >= 0 && var < 3) {
        posMinus[var] -= delta;
        posPlus[var] += delta;
    }
    cav.apply(posPlus, centroid, tPlus, Eplus, Bplus);
    cav.apply(posMinus, centroid, tMinus, Eminus, Bminus);
    dE = (Eplus - Eminus)/2/delta;
    dB = (Bplus - Bminus)/2/delta;
    if (verbose) {
        std::cerr << "Partial " << var << " " << delta << std::endl;
        std::cerr << "R: " << posPlus << " " << posMinus << std::endl;
        std::cerr << "t: " << tPlus << " " << tMinus << std::endl;
        std::cerr << "B: " << Bplus << " " << Bminus << std::endl;
        std::cerr << "E: " << Eplus << " " << Eminus << std::endl;
    }
}

// curl B = 1/c^2 dE/dt
Vector_t testMaxwell4(VariableRFCavityFringeField& cav, Vector_t pos, double t, double deltaPos, double deltaT) {
    bool verbose = false;
    Vector_t dummy;
    Vector_t dBdx(0, 0, 0);
    Vector_t dBdy(0, 0, 0);
    Vector_t dBdz(0, 0, 0);
    Vector_t dEdt(0, 0, 0);
    partial(cav, pos, t, deltaPos, 0, dummy, dBdx);
    partial(cav, pos, t, deltaPos, 1, dummy, dBdy);
    partial(cav, pos, t, deltaPos, 2, dummy, dBdz);
    partial(cav, pos, t, deltaT,   3, dEdt, dummy);
    Vector_t curlB(
        dBdy[2] - dBdz[1],
        dBdz[0] - dBdx[2],
        dBdx[1] - dBdy[0]
    );
    double c_l = Physics::c*1e-6; // 3e8 m/s = 300 mm/ns
    Vector_t result = dEdt - curlB*1e-1*c_l*c_l;

    if (verbose) {
        std::cerr << "dBdx           " << dBdx*1e-3 << std::endl;
        std::cerr << "dBdy           " << dBdy*1e-3 << std::endl;
        std::cerr << "dBdz           " << dBdz*1e-3 << std::endl;
        std::cerr << "dEdt           " << dEdt << std::endl;
        std::cerr << "c^2 curlB      " << curlB*1e-3*c_l*c_l << std::endl;
        std::cerr << "dEdt-c^2 curlB " << result << std::endl << std::endl;
    }
    return result;
}

Vector_t testMaxwell3(VariableRFCavityFringeField& cav, Vector_t pos, double t, double deltaPos, double deltaT) {
    bool verbose = false;
    Vector_t dummy;
    Vector_t dEdx(0, 0, 0);
    Vector_t dEdy(0, 0, 0);
    Vector_t dEdz(0, 0, 0);
    Vector_t dBdt(0, 0, 0);
    partial(cav, pos, t, deltaPos, 0, dEdx, dummy);
    partial(cav, pos, t, deltaPos, 1, dEdy, dummy);
    partial(cav, pos, t, deltaPos, 2, dEdz, dummy);
    partial(cav, pos, t, deltaT,   3, dummy, dBdt);
    Vector_t curlE(
        dEdy[2] - dEdz[1],
        dEdz[0] - dEdx[2],
        dEdx[1] - dEdy[0]
    );
    Vector_t result = dBdt*1e-1 + curlE;

    if (verbose) {
        std::cerr << "maxwell3 at R: " << pos << " t: " << t << " with dx: "
                  << deltaPos << " dt: " << deltaT << std::endl;
        std::cerr << "  dEdx           " << dEdx << std::endl;
        std::cerr << "  dEdy           " << dEdy << std::endl;
        std::cerr << "  dEdz           " << dEdz << std::endl;
        std::cerr << "  dBdt           " << dBdt << std::endl;
        std::cerr << "  curlE          " << curlE << std::endl;
        std::cerr << "  dBdt-curlE     " << result << std::endl << std::endl;
    }
    return result;
}


std::vector<double> testMaxwell1and2(VariableRFCavityFringeField& cav, Vector_t pos, double t, double deltaPos) {
    bool verbose = false;
    Vector_t dummy;
    Vector_t dEdx(0, 0, 0);
    Vector_t dEdy(0, 0, 0);
    Vector_t dEdz(0, 0, 0);
    Vector_t dBdx(0, 0, 0);
    Vector_t dBdy(0, 0, 0);
    Vector_t dBdz(0, 0, 0);
    partial(cav, pos, t, deltaPos, 0, dEdx, dBdx);
    partial(cav, pos, t, deltaPos, 1, dEdy, dBdy);
    partial(cav, pos, t, deltaPos, 2, dEdz, dBdz);

    double divE = dEdx[0] + dEdy[1] + dEdz[2];
    double divB = dBdx[0] + dBdy[1] + dBdz[2];

    if (verbose) {
        std::cerr << "maxwell1+2 at R: " << pos << " t: " << t << " with dx: "
                  << deltaPos << std::endl;
        std::cerr << "  dEidi          "
                  << dEdx[0] << " " << dEdy[1] << " " << dEdz[2] << std::endl;
        std::cerr << "  dBidi          "
                  << dBdx[0] << " " << dBdy[1] << " " << dBdz[2] << std::endl;
        std::cerr << "  divE           " << divE << std::endl;
        std::cerr << "  divB           " << divB << std::endl;
    }
    std::vector<double> result = {divE, divB};
    return result;
}


TEST_F(VariableRFCavityFringeFieldTest, TestField) {
    Vector_t centroid(0., 0., 0.);
    double t = 0.;
    cav2.setMaxOrder(4);
    std::cerr << "\nOff midplane, 45 degree phase, in fringe field" << std::endl;
    std::cerr << "order B        E   max1   max2   maxwell3    maxwell4" << std::endl;
    t = 0.125;
    for (double s = 0; s < 1000.; s += 10.) {
        Vector_t B0(0., 0., 0.);
        Vector_t E0(0., 0., 0.);
        Vector_t B10(0., 0., 0.);
        Vector_t E10(0., 0., 0.);
        Vector_t R(0., 0., 0.);
        R = Vector_t(0., 0., s);
        cav2.apply(R, centroid, t, E0, B0);
        R = Vector_t(0., 10., s);
        cav2.apply(R, centroid, t, E10, B10);
        std::cerr << s << " " << E0 << " " << E10 << " " << B10 << std::endl;
    }
}

TEST_F(VariableRFCavityFringeFieldTest, TestMaxwell) {
    //double pi = Physics::pi;
    Vector_t centroid(0., 0., 0.);
    Vector_t B(0., 0., 0.);
    Vector_t E(0., 0., 0.);
    Vector_t R(0., 0., 500.);
    double t = 0.;
    std::cerr << "\nOff midplane, 45 degree phase, in fringe field" << std::endl;
    std::cerr << "order B        E   max1   max2   maxwell3    maxwell4" << std::endl;
    R = Vector_t(0., 1., 750.);
    t = 0.125;
    for (size_t i = 0; i < 10; ++i) {
        cav2.setMaxOrder(i);
        cav2.apply(R, centroid, t, E, B);
        Vector_t result1 = testMaxwell3(cav2, R, t, 0.01, 0.0001);
        Vector_t result2 = testMaxwell4(cav2, R, t, 0.01, 0.0001);
        std::vector<double> div = testMaxwell1and2(cav2, R, t, 0.01);
        std::cerr << i << " ** " << B << " " << E << " " << div[0] << " "
                  << div[1] << " " << result1 << " " << result2 << std::endl;

        if (i > 0) {
            EXPECT_LT(div[0], 1e-3);
            EXPECT_LT(div[1], 1e-3);
            EXPECT_LT(euclidean_norm(result1), 1e-3);
            EXPECT_LT(euclidean_norm(result2), 1e-3);
        }
    }
}

TEST_F(VariableRFCavityFringeFieldTest, TestOrder) {
    Vector_t centroid(0., 0., 0.);
    Vector_t B(0., 0., 0.);
    Vector_t E(0., 0., 0.);
    Vector_t R(0., 0., 500.);
    double t = 0.;
    for (size_t i = 0; i < 20; ++i) {
        std::cerr << "Max Order " << i << std::endl;
        cav2.setMaxOrder(i);
        cav2.apply(R, centroid, t, E, B);
    }
}
