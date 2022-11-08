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
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/EndFieldModel/Enge.h"
#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"
#include "Elements/OpalAsymmetricEnge.h"

void setReal(OpalAsymmetricEnge& enge, std::string attName, double value) {
    Attribute* att = enge.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

void setRealArray(OpalAsymmetricEnge& enge,
                  std::string attName,
                  std::vector<double> value) {
    Attribute* att = enge.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setRealArray(*att, value);
}

TEST(OpalAsymmetricEngeTest, TestUpdate) {
    using namespace endfieldmodel;
    OpalAsymmetricEnge enge;
    enge.setOpalName("TEST_AENGE");
    setReal(enge, "X0_START", 1.0);
    setReal(enge, "LAMBDA_START", 2.0);
    setRealArray(enge, "COEFFICIENTS_START", {3.0, 4.0});
    setReal(enge, "X0_END", -1.0);
    setReal(enge, "LAMBDA_END", -2.0);
    setRealArray(enge, "COEFFICIENTS_END", {-3.0, -4.0, -5.0});
    enge.update();
    std::shared_ptr<EndFieldModel> myEFM =
                                EndFieldModel::getEndFieldModel("TEST_AENGE");
    AsymmetricEnge* myEnge = dynamic_cast<AsymmetricEnge*>(myEFM.get());
    std::shared_ptr<Enge> engeStart = myEnge->getEngeStart();
    std::shared_ptr<Enge> engeEnd = myEnge->getEngeEnd();
    EXPECT_EQ(myEnge->getX0Start(), 1.0);
    EXPECT_EQ(engeStart->getLambda(), 2.0);
    ASSERT_EQ(engeStart->getCoefficients().size(), 2u);
    EXPECT_EQ(engeStart->getCoefficients()[0], 3.0);
    EXPECT_EQ(engeStart->getCoefficients()[1], 4.0);
    EXPECT_EQ(myEnge->getX0End(), -1.0);
    EXPECT_EQ(engeEnd->getLambda(), -2.0);
    ASSERT_EQ(engeEnd->getCoefficients().size(), 3u);
    EXPECT_EQ(engeEnd->getCoefficients()[0], -3.0);
    EXPECT_EQ(engeEnd->getCoefficients()[1], -4.0);
    EXPECT_EQ(engeEnd->getCoefficients()[2], -5.0);
}