//
// Unit tests for class VariableRFCavity
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

#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Physics/Physics.h"

#include "gtest/gtest.h"

#include <vector>

#include "Physics/Units.h"

void testNull(VariableRFCavity& cav1) {
    std::shared_ptr<AbstractTimeDependence> null_poly(nullptr);
    EXPECT_DOUBLE_EQ(cav1.getLength(), 0.);
    EXPECT_EQ(cav1.getAmplitudeModel(), null_poly);
    EXPECT_EQ(cav1.getPhaseModel(), null_poly);
    EXPECT_EQ(cav1.getFrequencyModel(), null_poly);
}

TEST(VariableRFCavityTest, TestConstructorEtc) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    EXPECT_EQ(cav1.getName(), "");
    testNull(cav1);
    VariableRFCavity cav2("a_name");
    EXPECT_EQ(cav2.getName(), "a_name");
    testNull(cav1);
    // and now we implicitly check the destructor doesnt throw up on
    // case where everything is initialised to nullptr
}

// Is the obscure pointer-to-member function syntax appropriate here? I have
// been doing too much python where this stuff is easy
void testGetSet(VariableRFCavity& cav1,
                std::shared_ptr<AbstractTimeDependence> (VariableRFCavity::*getMethod)() const,
                void (VariableRFCavity::*setMethod)(std::shared_ptr<AbstractTimeDependence>)) {
    std::shared_ptr<AbstractTimeDependence> poly_1(new PolynomialTimeDependence(std::vector<double>(1, 1.)));
    std::shared_ptr<AbstractTimeDependence> poly_2(new PolynomialTimeDependence(std::vector<double>(2, 2.)));

    (cav1.*setMethod)(poly_1);
    EXPECT_EQ((cav1.*getMethod)(), poly_1);  // shallow equals is okay
    (cav1.*setMethod)(poly_2);
    EXPECT_EQ((cav1.*getMethod)(), poly_2);  // shallow equals is okay
    (cav1.*setMethod)(poly_2);
    EXPECT_EQ((cav1.*getMethod)(), poly_2);  // shallow equals is okay
    (cav1.*setMethod)(nullptr);  // and this deletes the memory
}

TEST(VariableRFCavityTest, TestGetSet) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    testGetSet(cav1,
               &VariableRFCavity::getAmplitudeModel,
               &VariableRFCavity::setAmplitudeModel);
    testGetSet(cav1,
               &VariableRFCavity::getPhaseModel,
               &VariableRFCavity::setPhaseModel);
    testGetSet(cav1,
               &VariableRFCavity::getFrequencyModel,
               &VariableRFCavity::setFrequencyModel);
    testNull(cav1);
    cav1.setLength(99.);
    EXPECT_DOUBLE_EQ(cav1.getLength(), 99.);
}

TEST(VariableRFCavityTest, TestAssignmentNull) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    VariableRFCavity cav2;
    cav2.getLength();  // stop compiler "optimising" to copy constructor
    cav2 = cav1;
    testNull(cav2);
    VariableRFCavity cav3(cav2);
    testNull(cav3); // now this is really the copy constructor
}

TEST(VariableRFCavityTest, TestAssignmentValue) {
    OpalTestUtilities::SilenceTest silencer;

    std::shared_ptr<AbstractTimeDependence> poly1(new PolynomialTimeDependence(std::vector<double>(1, 1.)));
    std::shared_ptr<AbstractTimeDependence> poly2(new PolynomialTimeDependence(std::vector<double>(1, 2.)));
    std::shared_ptr<AbstractTimeDependence> poly3(new PolynomialTimeDependence(std::vector<double>(1, 3.)));
    VariableRFCavity cav1;
    cav1.setPhaseModel(poly1);
    cav1.setAmplitudeModel(poly2);
    cav1.setFrequencyModel(poly3);
    cav1.setLength(99.);
    VariableRFCavity cav2(cav1);
    EXPECT_EQ(cav1.getPhaseModel()->getValue(1.),
              cav2.getPhaseModel()->getValue(1.));
    EXPECT_EQ(cav1.getAmplitudeModel()->getValue(1.),
              cav2.getAmplitudeModel()->getValue(1.));
    EXPECT_EQ(cav1.getFrequencyModel()->getValue(1.),
              cav2.getFrequencyModel()->getValue(1.));
    EXPECT_DOUBLE_EQ(cav1.getLength(), cav2.getLength());
}

TEST(VariableRFCavityTest, TestClone) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    cav1.setLength(99.);
    VariableRFCavity* cav2 = dynamic_cast<VariableRFCavity*>(cav1.clone());
    EXPECT_DOUBLE_EQ(cav1.getLength(), cav2->getLength());
    delete cav2;
}

TEST(VariableRFCavityTest, TestInitialiseFinalise) {
    OpalTestUtilities::SilenceTest silencer;

    // nothing to do here
}

TEST(VariableRFCavityTest, TestGetGeometry) {
    OpalTestUtilities::SilenceTest silencer;
    VariableRFCavity cav1;
    const VariableRFCavity& cav2(cav1);
    EXPECT_EQ(&cav1.getGeometry(), &cav2.getGeometry());
    cav1.setLength(99.);
    EXPECT_EQ(cav1.getGeometry().getElementLength(), cav1.getLength());
}

TEST(VariableRFCavityTest, TestBends) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    EXPECT_FALSE(cav1.bends());
}

TEST(VariableRFCavityTest, TestApplyField) {
    OpalTestUtilities::SilenceTest silencer;

    VariableRFCavity cav1;
    std::vector<double>  vec1;
    vec1.push_back(1.);
    vec1.push_back(2.);
    std::vector<double>  vec2;
    vec2.push_back(3.);
    vec2.push_back(4.);
    std::vector<double>  vec3;
    vec3.push_back(5.);
    vec3.push_back(6.);
    std::shared_ptr<AbstractTimeDependence> poly1(new PolynomialTimeDependence(vec1));
    std::shared_ptr<AbstractTimeDependence> poly2(new PolynomialTimeDependence(vec2));
    std::shared_ptr<AbstractTimeDependence> poly3(new PolynomialTimeDependence(vec3));
    cav1.setAmplitudeModel(poly1);
    cav1.setFrequencyModel(poly2);
    cav1.setPhaseModel(poly3);
    cav1.setLength(2.);
    cav1.setWidth(3.);
    cav1.setHeight(4.);
    Vector_t R({1., 1., 1.});
    Vector_t centroid({0., 0., 0.});
    Vector_t B({0., 0., 0.});
    Vector_t E({0., 0., 0.});
    for (double t = 0.0; t < 10.0e-9; t += 1.0e-9) {
        double phase = poly3->getValue(t);
        double amplitude = poly1->getValue(t);
        double integralF = poly2->getIntegral(t) * Units::MHz2Hz;
        double e_test = amplitude * sin(Physics::two_pi * integralF + phase);
        ASSERT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
        EXPECT_NEAR(0., E[0], 1.e-6);
        EXPECT_NEAR(0., E[1], 1.e-6);
        EXPECT_NEAR(e_test, E[2], 1.e-6);
        EXPECT_NEAR(0., B[0], 1.e-6);
        EXPECT_NEAR(0., B[1], 1.e-6);
        EXPECT_NEAR(0., B[2], 1.e-6);
    }
}

TEST(VariableRFCavityTest, TestApplyBoundingBox) {
    OpalTestUtilities::SilenceTest silencer;
    VariableRFCavity cav1;
    std::shared_ptr<AbstractTimeDependence> poly1(new PolynomialTimeDependence(std::vector<double>(1, 1.)));
    std::shared_ptr<AbstractTimeDependence> poly2(new PolynomialTimeDependence(std::vector<double>(2, 2.)));
    std::shared_ptr<AbstractTimeDependence> poly3(new PolynomialTimeDependence(std::vector<double>(3, 3.)));
    cav1.setAmplitudeModel(poly1);
    cav1.setFrequencyModel(poly2);
    cav1.setPhaseModel(poly3);
    cav1.setLength(2.);
    cav1.setHeight(3.);
    cav1.setWidth(4.);
    Vector_t R({0., 0., 1.});
    Vector_t centroid({0., 0., 0.});
    Vector_t B({0., 0., 0.});
    Vector_t E({0., 0., 0.});
    double t = 0;
    EXPECT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[2] = 2.-1e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[2] = 1.e-9;
    EXPECT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[2] = -1.e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[2] = 2.+1.e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[2] = 1.;
    R[1] = -1.5-1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[1] = +1.5+1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[1] = 0.;
    EXPECT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[0] = -2.-1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[0] = +2.+1e-9;
    EXPECT_TRUE(cav1.apply(R, Vector_t(0.0), t, E, B));
    R[0] = 0.;
    EXPECT_FALSE(cav1.apply(R, Vector_t(0.0), t, E, B));
}
