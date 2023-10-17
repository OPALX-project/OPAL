//
// Unit tests for OpalOffset element definition
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
#include "AbsBeamline/Offset.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Object.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalOffset/OpalGlobalCartesianOffset.h"
#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"

#include "gtest/gtest.h"

#include <cmath>

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
    Vector_t refDir(0.1, std::sqrt(0.99), 0.0);
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
    EXPECT_NEAR(pos[0], refPos[0], 1e-9);
    EXPECT_NEAR(pos[1], refPos[1], 1e-9);

    Vector_t dir = off->getEndDirection();
    EXPECT_NEAR(dir[0], refDir[0], 1e-9);
    EXPECT_NEAR(dir[1], refDir[1], 1e-9);

    EXPECT_FALSE(off->getIsLocal());
}

TEST(OpalOffsetTest, TestLocalCartesianOffset) {
    OpalLocalCartesianOffset offset;
    Vector_t refDir(0.1, std::sqrt(0.99), 0.0);
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
    EXPECT_NEAR(pos[0], refPos[0], 1e-9);
    EXPECT_NEAR(pos[1], refPos[1], 1e-9);

    Vector_t dir = off->getEndDirection();
    EXPECT_NEAR(dir[0], refDir[0], 1e-9);
    EXPECT_NEAR(dir[1], refDir[1], 1e-9);

    EXPECT_TRUE(off->getIsLocal());
}
