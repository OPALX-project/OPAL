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
#include "Elements/OpalEnge.h"

void setReal(OpalEnge& enge, std::string attName, double value) {
    Attribute* att = enge.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

void setRealArray(OpalEnge& enge,
                  std::string attName,
                  std::vector<double> value) {
    Attribute* att = enge.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setRealArray(*att, value);
}


TEST(OpalEngeTest, TestUpdate) {
    using namespace endfieldmodel;
    OpalEnge enge;
    enge.setOpalName("TEST_ENGE");
    setReal(enge, "X0", 1.0);
    setReal(enge, "LAMBDA", 2.0);
    setRealArray(enge, "COEFFICIENTS", {3.0, 4.0});
    enge.update();
    std::shared_ptr<EndFieldModel> myEFM =
                                EndFieldModel::getEndFieldModel("TEST_ENGE");
    Enge* myEnge = dynamic_cast<Enge*>(myEFM.get());
    EXPECT_EQ(myEnge->getX0(), 1.0);
    EXPECT_EQ(myEnge->getLambda(), 2.0);
    ASSERT_EQ(myEnge->getCoefficients().size(), 2u);
    EXPECT_EQ(myEnge->getCoefficients()[0], 3.0);
    EXPECT_EQ(myEnge->getCoefficients()[1], 4.0);
}