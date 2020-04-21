//
// Test PythonExprTest
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
#include <set>
#include <string>

#include "Util/Types.h"
#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"
#include "Expression/PythonExpr.h"

#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"


namespace {

    // The fixture for testing class Foo.
    class PythonExprTest : public ::testing::Test {
    protected:

        PythonExprTest() {
            // You can do set-up work for each test here.
        }

        virtual ~PythonExprTest() {
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


    TEST_F(PythonExprTest, EvaluatePythonExpr) {

        variableDictionary_t vars;
        double expected = 1.0;
        expected *= 2.0;

        functionDictionary_t funcs;
        client::function::type python;
        python = PythonExpression();
        funcs.insert(std::pair<std::string, client::function::type>
                ("python", python));

        std::string testexpr = "python(\"resources/test.py\", 1.0)";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));
        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        ASSERT_EQ(expected, boost::get<0>(result));
        ASSERT_TRUE(boost::get<1>(result));
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
