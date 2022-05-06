/*
 *  Copyright (c) 2015, Chris Rogers
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

#include "Fields/Interpolation/SolveFactory.h"
#include "Fields/Interpolation/SquarePolynomialVector.h"
#include "Fields/Interpolation/MMatrix.h"

#include "opal_test_utilities/SilenceTest.h"

using namespace interpolation;

class SolveFactoryTestFixture : public ::testing::Test {
  protected:


};

TEST(SolveFactoryTest, TestSolveNoDerivs) {
    OpalTestUtilities::SilenceTest silencer;

    // we make a reference poly vector
    std::vector<double> data(27);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;
    MMatrix<double> refCoeffs(3, 9, &data[0]);
    SquarePolynomialVector ref(2, refCoeffs);
    // Make a set of points
    std::vector< std::vector<double> > positions;
    std::vector< std::vector<double> > values;
    std::vector<double> pos(2);
    std::vector<double> val(3);
    for (pos[0] = 0.; pos[0] < 2.5; pos[0] += 1.)
        for (pos[1] = 0.; pos[1] < 2.5; pos[1] += 1.) {
            ref.F(&pos[0], &val[0]);
            positions.push_back(pos);
            values.push_back(val);
        }
    // now try to solve for the test poly vector; should be that the test and
    // ref are identical
    std::vector<std::vector<double> > deriv_pos;
    std::vector< std::vector<int> > deriv_index;
    SolveFactory fac(2, 2, 3, positions, deriv_pos, deriv_index);
    SquarePolynomialVector* vec = fac.PolynomialSolve(values,
                                                      deriv_pos);
    MMatrix<double> testCoeffs = vec->GetCoefficientsAsMatrix();
    ASSERT_EQ(testCoeffs.num_row(), (size_t)3);
    ASSERT_EQ(testCoeffs.num_col(), (size_t)9);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6);
        }
    }
    delete vec;
}

TEST(SolveFactoryTest, TestSolveDerivs) {
    // we make a reference poly vector
    std::vector<double> data(27);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;
    MMatrix<double> refCoeffs(3, 9, &data[0]);
    SquarePolynomialVector ref(2, refCoeffs);
    std::vector<int> derivIndex01({0, 1});
    SquarePolynomialVector refDeriv01 = ref.Deriv(&derivIndex01[0]);
    std::vector<int> derivIndex11({1, 1});
    SquarePolynomialVector refDeriv11 = ref.Deriv(&derivIndex11[0]);
    std::vector<int> derivIndex10({1, 0});
    SquarePolynomialVector refDeriv10 = ref.Deriv(&derivIndex10[0]);
    // Make a set of points
    std::vector< std::vector<double> > positions;
    std::vector< std::vector<double> > values;
    std::vector<double> pos(2);
    std::vector<double> val(3);
    for (pos[0] = 0.; pos[0] < 1.5; pos[0] += 1.) {
        for (pos[1] = 0.; pos[1] < 1.5; pos[1] += 1.) {
            ref.F(&pos[0], &val[0]);
            positions.push_back(pos);
            values.push_back(val);
        }
    }

    // fill by hand the
    std::vector< std::vector<double> > deriv_pos( {
        {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}
    });
    std::vector< std::vector<int> > deriv_indices( {
        {0,1}, {0,1}, {1,1}, {1,0}, {1,0}
    });
    std::vector<std::vector<double> > deriv_values(deriv_pos.size(),
                                                   std::vector<double>(3));
    for (size_t i = 0; i < deriv_pos.size(); ++i) {
        SquarePolynomialVector refDeriv = ref.Deriv(&deriv_indices[i][0]);
        refDeriv.F(&deriv_pos[i][0], &deriv_values[i][0]);
    }

    // now try to solve for the test poly vector; should be that the test and
    // ref are identical
    SquarePolynomialVector* vec = nullptr;
    try {
        SolveFactory fac(2, 2, 3, positions, deriv_pos, deriv_indices);
        vec = fac.PolynomialSolve(values, deriv_values);
    } catch (ClassicException& exc) {
        std::cerr << exc.what() << std::endl;
        ASSERT_TRUE(false) << "Should not throw";
    }
    MMatrix<double> testCoeffs = vec->GetCoefficientsAsMatrix();
    ASSERT_EQ(testCoeffs.num_row(), (size_t)3);
    ASSERT_EQ(testCoeffs.num_col(), (size_t)9);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6);
}
