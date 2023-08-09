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

#include "Physics/Units.h"
#include "Attributes/Attributes.h"
#include "AbsBeamline/Offset.h"
#include "Elements/OpalOffset/OpalGlobalCartesianOffset.h"
#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"

using OpalOffset::OpalGlobalCartesianOffset;
using OpalOffset::OpalLocalCartesianOffset;

void setReal(OpalGlobalCartesianOffset& mag, std::string attName, double value) {
    Attribute* att = mag.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

void setReal(OpalLocalCartesianOffset& mag, std::string attName, double value) {
    Attribute* att = mag.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

TEST(OpalOffsetTest, TestGlobalCartesianOffset) {
    OpalGlobalCartesianOffset offset;
    Vector_t refDir(0.1, sqrt(0.99), 0.0);
    Vector_t refPos(1, 2, 0.0);
    setReal(offset, "END_POSITION_X", refPos[0]);
    setReal(offset, "END_POSITION_Y", refPos[1]);
    setReal(offset, "END_NORMAL_X", refDir[0]);
    setReal(offset, "END_NORMAL_Y", refDir[1]);
    offset.update();
    ElementBase* element = offset.getElement();
    Offset* off = dynamic_cast<Offset*>(element);
    ASSERT_NE(off, nullptr);

    Vector_t pos = off->getEndPosition();
    EXPECT_NEAR(pos[0], refPos[0]*Units::m2mm, 1e-9);
    EXPECT_NEAR(pos[1], refPos[1]*Units::m2mm, 1e-9);

    Vector_t dir = off->getEndDirection();
    EXPECT_NEAR(dir[0], refDir[0], 1e-9);
    EXPECT_NEAR(dir[1], refDir[1], 1e-9);

    EXPECT_FALSE(off->getIsLocal());
}

TEST(OpalOffsetTest, TestLocalCartesianOffset) {
    OpalLocalCartesianOffset offset;
    Vector_t refDir(0.1, sqrt(0.99), 0.0);
    Vector_t refPos(1, 2, 0.0);
    setReal(offset, "END_POSITION_X", refPos[0]);
    setReal(offset, "END_POSITION_Y", refPos[1]);
    setReal(offset, "END_NORMAL_X", refDir[0]);
    setReal(offset, "END_NORMAL_Y", refDir[1]);
    offset.update();
    ElementBase* element = offset.getElement();
    Offset* off = dynamic_cast<Offset*>(element);
    ASSERT_NE(off, nullptr);

    Vector_t pos = off->getEndPosition();
    EXPECT_NEAR(pos[0], refPos[0]*Units::m2mm, 1e-9);
    EXPECT_NEAR(pos[1], refPos[1]*Units::m2mm, 1e-9);

    Vector_t dir = off->getEndDirection();
    EXPECT_NEAR(dir[0], refDir[0], 1e-9);
    EXPECT_NEAR(dir[1], refDir[1], 1e-9);

    EXPECT_TRUE(off->getIsLocal());
}
