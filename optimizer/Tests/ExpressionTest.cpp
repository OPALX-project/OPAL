//
// Test ExpressionTest
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

#include "Expression/Expression.h"
#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"

#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"


namespace {


struct my_abs {

    Expressions::Result_t operator()(client::function::arguments_t args) {

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(fabs(value), true);
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
    class ExpressionTest : public ::testing::Test {
    protected:

        ExpressionTest() {
            // You can do set-up work for each test here.
        }

        virtual ~ExpressionTest() {
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

    TEST_F(ExpressionTest, ParseExpression) {

        EXPECT_NO_THROW({
            Expression("(1 + A * A * 2 * b + 2)");
        });

    }

    TEST_F(ExpressionTest, RequestedVars) {

        std::string testexpr = "abs(A * A * 2.0 * b + 2.0)";
        boost::scoped_ptr<Expression> e(new Expression(testexpr));

        std::set<std::string> reqVars = e->getReqVars();

        ASSERT_EQ(static_cast<size_t>(2), reqVars.size());

    }

    TEST_F(ExpressionTest, EvaluateExpression) {

        functionDictionary_t funcs;
        client::function::type abs_func;
        abs_func = my_abs();
        funcs.insert(std::pair<std::string, client::function::type>
                ("abs", abs_func));

        std::string testexpr = "abs(1.0 + A * A * 2.0 * B + 2.0)";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));

        double a = 5.2;
        double b = -10.2;

        variableDictionary_t vars;
        vars.insert(std::pair<std::string, double>("A", a));
        vars.insert(std::pair<std::string, double>("B", b));

        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        double expected = fabs(1 + a * a * 2 * b + 2);
        ASSERT_DOUBLE_EQ(expected, boost::get<0>(result));
        ASSERT_TRUE(boost::get<1>(result));

    }

    TEST_F(ExpressionTest, EvaluateCombinedExpression) {

        functionDictionary_t funcs;
        client::function::type abs_func;
        abs_func = my_abs();
        funcs.insert(std::pair<std::string, client::function::type>
                ("abs", abs_func));

        std::string testexpr = "abs(1.0 + A * 2.0) + abs(B + 2.0)";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));

        double a = 5.2;
        double b = -10.2;

        variableDictionary_t vars;
        vars.insert(std::pair<std::string, double>("A", a));
        vars.insert(std::pair<std::string, double>("B", b));

        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        double expected = fabs(1 + a * 2) + fabs(b + 2);
        ASSERT_DOUBLE_EQ(expected, boost::get<0>(result));
        ASSERT_TRUE(boost::get<1>(result));

    }

    TEST_F(ExpressionTest, EvaluateNestedExpression) {

        functionDictionary_t funcs;
        client::function::type abs_func;
        abs_func = my_abs();
        funcs.insert(std::pair<std::string, client::function::type>
                ("abs", abs_func));
        client::function::type pow_func;
        pow_func = my_pow();
        funcs.insert(std::pair<std::string, client::function::type>
                ("pow", pow_func));

        std::string testexpr = "abs(1.0 + pow(A, 3) * 2.0) + abs(B + 2.0)";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));

        double a = 5.2;
        double b = -10.2;

        variableDictionary_t vars;
        vars.insert(std::pair<std::string, double>("A", a));
        vars.insert(std::pair<std::string, double>("B", b));

        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        double expected = fabs(1 + pow(a, 3.0) * 2) + fabs(b + 2);
        ASSERT_DOUBLE_EQ(expected, boost::get<0>(result));
        ASSERT_TRUE(boost::get<1>(result));

    }

    TEST_F(ExpressionTest, EvaluateBooleanExpression) {

        std::string testexpr = "a > 5.2 - 1e-6";
        boost::scoped_ptr<Expression> e(new Expression(testexpr));

        double a = 5.2;

        variableDictionary_t vars;
        vars.insert(std::pair<std::string, double>("a", a));

        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });

        bool expected = true;
        ASSERT_DOUBLE_EQ(expected, boost::get<0>(result));
        ASSERT_TRUE(boost::get<1>(result));

        Expressions::OperatorType_t op = e->getOpType();
        ASSERT_EQ(Expressions::INEQ_RHS, op);
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}