//
// Unit tests for class Undulator
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#include <vector>

#include "gtest/gtest.h"

#include "AbsBeamline/Undulator.h"
#include "BeamlineCore/UndulatorRep.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/Vektor.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"

void testNull(Undulator& und) {

    EXPECT_EQ(und.getTypeString(), "Undulator");
    EXPECT_DOUBLE_EQ(und.getK(), 0.);
    EXPECT_DOUBLE_EQ(und.getLambda(), 0.);
    EXPECT_EQ(und.getNumPeriods(), (unsigned int) 0);
    EXPECT_DOUBLE_EQ(und.getAngle(), 0.);
    EXPECT_EQ(und.getFilename(), "");
    std::vector<double> nullVec (3, 0.0);
    std::vector<double> meshLength =  und.getMeshLength();
    std::vector<double> meshResolution =  und.getMeshResolution();
    EXPECT_EQ(meshLength, nullVec);
    EXPECT_EQ(meshResolution, nullVec);
    EXPECT_EQ(und.getTruncationOrder(), (unsigned int) 2);
    EXPECT_DOUBLE_EQ(und.getTotalTime(), 0.);
    EXPECT_DOUBLE_EQ(und.getDtBunch(), 0.);
    EXPECT_EQ(und.getHasBeenSimulated(), false);
}

TEST(UndulatorTest, TestConstructorAndGets) {
    OpalTestUtilities::SilenceTest silencer;

    UndulatorRep und1;
    EXPECT_EQ(und1.getName(), "");
    testNull(und1);
    UndulatorRep und2("a_name");
    EXPECT_EQ(und2.getName(), "a_name");
    testNull(und2);
    UndulatorRep und3(und2);
    EXPECT_EQ(und3.getName(), und2.getName());
    testNull(und3);
}

TEST(UndulatorTest, TestBends) {
    OpalTestUtilities::SilenceTest silencer;

    UndulatorRep und1;
    EXPECT_FALSE(und1.bends());
}

TEST(UndulatorTest, TestGetSet) {
    OpalTestUtilities::SilenceTest silencer;

    UndulatorRep und1;
    
    und1.setK(1.0);
    EXPECT_DOUBLE_EQ(und1.getK(), 1.0);
    und1.setLambda(1.0);
    EXPECT_DOUBLE_EQ(und1.getLambda(), 1.0);
    und1.setNumPeriods(1);
    EXPECT_EQ(und1.getNumPeriods(), (unsigned int) 1);
    und1.setAngle(45.0);
    EXPECT_DOUBLE_EQ(und1.getAngle(), 45.0);
    und1.setFilename("file");
    EXPECT_EQ(und1.getFilename(), "file");
    std::vector<double> testVec (3, 1.0);
    und1.setMeshLength(testVec);
    EXPECT_EQ(und1.getMeshLength(), testVec);
    und1.setMeshResolution(testVec);
    EXPECT_EQ(und1.getMeshResolution(), testVec);
    und1.setTruncationOrder(1);
    EXPECT_EQ(und1.getTruncationOrder(), (unsigned int) 1);
    und1.setTotalTime(1.0);
    EXPECT_DOUBLE_EQ(und1.getTotalTime(), 1.0);
    und1.setDtBunch(1.0);
    EXPECT_DOUBLE_EQ(und1.getDtBunch(), 1.0);
    und1.setHasBeenSimulated(true);
    EXPECT_EQ(und1.getHasBeenSimulated(), true);
}

TEST(UndulatorTest, TestApplyFullWaveSolver) {
    OpalTestUtilities::SilenceTest silencer;

    UndulatorRep und;
    PartData* partData = new PartData();
    PartBunch* bunch = new PartBunch(partData);
    CoordinateSystemTrafo refToLocalCSTrafo;
    
    und.apply(bunch, refToLocalCSTrafo);
    EXPECT_EQ(und.getHasBeenSimulated(), true);
}