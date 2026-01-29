//
// Test HashNameGeneratorTest
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
#include "Util/NativeHashGenerator.h"
#include "gtest/gtest.h"
#include <fstream>

namespace {

    // The fixture for testing class Foo.
    class HashNameGeneratorTest : public ::testing::Test {
    protected:

        HashNameGeneratorTest() {
            // You can do set-up work for each test here.
        }

        virtual ~HashNameGeneratorTest() {
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
    };

    TEST_F(HashNameGeneratorTest, HashName) {

        std::vector<std::string> params;
        params.push_back("transferline");
        params.push_back("sigmax=5.05");
        params.push_back("sigmay=6.05");

        std::string hash = NativeHashGenerator::generate(params);
        std::string expected = "88a02dc533c53bc156f54a4fcc54397b";
        ASSERT_STREQ(expected.c_str(), hash.c_str());
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}