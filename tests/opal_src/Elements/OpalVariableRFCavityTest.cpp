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

#include "Attributes/Attributes.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Elements/OpalVariableRFCavity.h"
#include "BeamlineCore/DriftRep.h"

#include "opal_test_utilities/SilenceTest.h"

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
