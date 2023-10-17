//
// Unit tests for OpalVariableRFCavity element definition
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
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "Elements/OpalVariableRFCavity.h"

#include "gtest/gtest.h"

TEST(OpalVariableRFCavityTest, TestConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    OpalVariableRFCavity cav1;
    EXPECT_EQ((&cav1)->getOpalName(), "VARIABLE_RF_CAVITY");
    OpalVariableRFCavity cav2("name", &cav1);
    EXPECT_EQ((&cav2)->getOpalName(), "name");
    OpalVariableRFCavity* cav3 = cav2.clone();
    EXPECT_EQ(cav3->getOpalName(), "name");
    OpalVariableRFCavity* cav4 = cav2.clone("other_name");
    EXPECT_EQ(cav4->getOpalName(), "other_name");

    delete cav4;
    delete cav3;
}

TEST(OpalVariableRFCavityTest, TestFillRegisteredAttributes) {
    OpalTestUtilities::SilenceTest silencer;

/*
    PolynomialTimeDependence* pd1 = new PolynomialTimeDependence();
    PolynomialTimeDependence* pd2 = new PolynomialTimeDependence();
    PolynomialTimeDependence* pd3 = new PolynomialTimeDependence();
    AbstractTimeDependence::setTimeDependence("pd1", pd1);
    AbstractTimeDependence::setTimeDependence("pd2", pd2);
    AbstractTimeDependence::setTimeDependence("pd3", pd3);

    VariableRFCavity cav("my_name");
    cav.setLength(99.);
    cav.setPhaseModel(pd1);
    cav.setAmplitudeModel(pd2);
    cav.setFrequencyModel(pd3);

    OpalVariableRFCavity opal_cav;
    OpalVariableRFCavity parent;  // dummy parent to prevent segv
    opal_cav.setParent(&parent);
    opal_cav.fillRegisteredAttributes(cav);

    Attribute* null_att = nullptr;
    EXPECT_EQ(opal_cav.findAttribute("NONSENSE ATTRIBUTE ASDASDA"), null_att);
    ASSERT_NE(opal_cav.findAttribute("L"), null_att);
    EXPECT_EQ(Attributes::getReal(*opal_cav.findAttribute("L")), 99.);
    ASSERT_NE(opal_cav.findAttribute("PHASE_MODEL"), null_att);
    EXPECT_EQ(Attributes::getString(*opal_cav.findAttribute("PHASE_MODEL")),
              "pd1");
    ASSERT_NE(opal_cav.findAttribute("AMPLITUDE_MODEL"), null_att);
    EXPECT_EQ(Attributes::getString(*opal_cav.findAttribute("AMPLITUDE_MODEL")),
              "pd2");
    ASSERT_NE(opal_cav.findAttribute("FREQUENCY_MODEL"), null_att);
    EXPECT_EQ(Attributes::getString(*opal_cav.findAttribute("FREQUENCY_MODEL")),
              "pd3");

    // try to fill a VariableRFCavity using an ElementBase that is not a
    // VariableRFCavity - should throw
    DriftRep drift("test");
    EXPECT_THROW(opal_cav.fillRegisteredAttributes(drift),
                 OpalException);
*/
}

TEST(OpalVariableRFCavityTest, TestUpdate) {
    OpalTestUtilities::SilenceTest silencer;

//    EXPECT_TRUE(false);
}
