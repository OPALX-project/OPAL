//
// Test CmdArgumentsTest
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#include "Util/CmdArguments.h"
#include "gtest/gtest.h"
#include "boost/smart_ptr.hpp"

namespace {

    // The fixture for testing class Foo.
    class CmdArgumentsTest : public ::testing::Test {
    protected:

        CmdArgumentsTest() {
            // You can do set-up work for each test here.
            int argc        = 3;
            char exe_name[] = "test";
            char a1[]       = "--arg1=val1";
            char a2[]       = "--arg2=2.2";
            char *argv[]    = { exe_name, a1, a2 };

            args_.reset(new CmdArguments(argc, argv));
        }

        virtual ~CmdArgumentsTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case
        boost::scoped_ptr<CmdArguments> args_;
    };

    TEST_F(CmdArgumentsTest, RetrieveCorrectFatal) {

        std::string arg1 = args_->getArg<std::string>("arg1", true);
        double arg2 = args_->getArg<double>("arg2", true);

        EXPECT_EQ("val1", arg1) << "first argument string value wrong";
        EXPECT_EQ(2.2, arg2)    << "second argument double value wrong";
    }

    TEST_F(CmdArgumentsTest, ThrowOnNotPresentFatal) {

        EXPECT_ANY_THROW(
            args_->getArg<std::string>("arg11", true)
        );
    }

    TEST_F(CmdArgumentsTest, CorrectDefaultIfNotPresent) {

        double arg = args_->getArg<double>("arg22", 10.0, false);
        EXPECT_EQ(10.0, arg) << "second argument double value wrong";
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
