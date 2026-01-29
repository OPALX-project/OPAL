//
// Test SumErrSqExpressionTest
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
#include "Expression/SumErrSq.h"

#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"


namespace {

    // The fixture for testing class Foo.
    class SumErrSqExpressionTest : public ::testing::Test {
    protected:

        SumErrSqExpressionTest() {
            // You can do set-up work for each test here.
        }

        virtual ~SumErrSqExpressionTest() {
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


    TEST_F(SumErrSqExpressionTest, EvaluateSumErrSqExpression) {

        variableDictionary_t vars;
        double expected = (3.087242557177229e-04*3.087242557177229e-04 +
                           3.127445619624299e-04*3.127445619624299e-04 +
                           3.185324887508158e-04*3.185324887508158e-04) / 3.0;
        expected = sqrt(expected);

        functionDictionary_t funcs;
        client::function::type errsumsq;
        errsumsq = SumErrSq();
        funcs.insert(std::pair<std::string, client::function::type>
                ("sumErrSq", errsumsq));

        std::string testexpr = "sumErrSq(\"resources/measurement_test.dat\", \"rms_x\", \"resources/test.stat\")";
        boost::scoped_ptr<Expression> e(new Expression(testexpr, funcs));
        Expressions::Result_t result;
        EXPECT_NO_THROW({
            result = e->evaluate(vars);
        });


        //XXX: expected uses the nearest (and NOT interpolated) rms_x values
        ASSERT_NEAR(expected, boost::get<0>(result), 1e-6);
        ASSERT_TRUE(boost::get<1>(result));
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
