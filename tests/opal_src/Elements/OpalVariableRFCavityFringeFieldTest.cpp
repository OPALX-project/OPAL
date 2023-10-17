//
// Unit tests for OpalVariableRFCavityFringeField element definition
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

#include "Elements/OpalVariableRFCavityFringeField.h"

#include "gtest/gtest.h"

TEST(OpalVariableRFCavityFringeFieldTest, TestConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    OpalVariableRFCavityFringeField cav1;
    EXPECT_EQ((&cav1)->getOpalName(), "VARIABLE_RF_CAVITY_FRINGE_FIELD");
    OpalVariableRFCavityFringeField cav2("name", &cav1);
    EXPECT_EQ((&cav2)->getOpalName(), "name");
    OpalVariableRFCavityFringeField* cav3 = cav2.clone();
    EXPECT_EQ(cav3->getOpalName(), "name");
    OpalVariableRFCavityFringeField* cav4 = cav2.clone("other_name");
    EXPECT_EQ(cav4->getOpalName(), "other_name");

    delete cav4;
    delete cav3;
}
