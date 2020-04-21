//
// Test FromFileExpressionTest
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

#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"

#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"


namespace {

struct my_sqrt {

    Expressions::Result_t operator()(client::function::arguments_t args) {

        double value = boost::get<double>(args[0]);
        if (value < 0.0) return boost::make_tuple(0.0, false);
        return boost::make_tuple(sqrt(value), true);
    }
};

struct my_pow {

    Expressions::Result_t operator()(client::function::arguments_t args) {

        double base = boost::get<double>(args[0]);
        double exponent = boost::get<double>(args[1]);
        return boost::make_tuple(pow(base, exponent), true);
    }
};

    // The fixture for testing class Foo.
    class FromFileExpressionTest : public ::testing::Test {
    protected:

        FromFileExpressionTest() {
            // You can do set-up work for each test here.
        }

        virtual ~FromFileExpressionTest() {
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


    TEST_F(FromFileExpressionTest, EvaluateFromFileExpression) {

        variableDictionary_t vars;
        double expected = 0.3126 + 0.3561 + 0.4242;

        functionDictionary_t funcs;
        client::function::type ff;
        ff = FromFile();
        funcs.insert(std::pair<std::string, client::function::type>
                ("fromFile", ff));

        std::string testexpr = "fromFile(\"resources/fromfile_test.dat\")";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));
        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        ASSERT_NEAR(expected, boost::get<0>(result), 1e-6);
        ASSERT_TRUE(boost::get<1>(result));
    }

    TEST_F(FromFileExpressionTest, EvaluateCombinedFromFileExpression) {

        variableDictionary_t vars;
        double expected = sqrt(pow(0.3126 + 0.3561 + 0.4242, 2) + pow(0.1263 - 0.5613 + 0.2424, 2));

        functionDictionary_t funcs;
        client::function::type ff;
        ff = FromFile();
        funcs.insert(std::pair<std::string, client::function::type>
                ("fromFile", ff));
        client::function::type sqrt_func;
        sqrt_func = my_sqrt();
        funcs.insert(std::pair<std::string, client::function::type>
                ("sqrt", sqrt_func));
        client::function::type pow_func;
        pow_func = my_pow();
        funcs.insert(std::pair<std::string, client::function::type>
                ("pow", pow_func));

        std::string testexpr = "sqrt(pow(fromFile(\"resources/fromfile_test.dat\"), 2) + "
            "pow(fromFile(\"resources/fromfile_test2.dat\"), 2))";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));
        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        ASSERT_NEAR(expected, boost::get<0>(result), 1e-6);
        ASSERT_TRUE(boost::get<1>(result));
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}