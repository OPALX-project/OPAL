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

#include "Fields/Interpolation/SquarePolynomialVector.h"
#include "Fields/Interpolation/PolynomialCoefficient.h"

#include "opal_test_utilities/SilenceTest.h"

using namespace interpolation;

TEST(SquarePolynomialVectorTest, TestConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    SquarePolynomialVector test1;
    EXPECT_EQ(test1.PointDimension(), (unsigned int)0);
    EXPECT_EQ(test1.ValueDimension(), (unsigned int)0);
    EXPECT_EQ(test1.GetCoefficientsAsMatrix().num_col(), (unsigned int)0);
    EXPECT_EQ(test1.GetCoefficientsAsMatrix().num_row(), (unsigned int)0);

    std::vector<double> data(18);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;
    MMatrix<double> refCoeffs(2, 9, &data[0]); // 2x9 -> c, x, y, xy, xx, xxy, xxyy, xyy, yy,
    SquarePolynomialVector test2(2, refCoeffs);
    MMatrix<double> testCoeffs = test2.GetCoefficientsAsMatrix();
    ASSERT_EQ(testCoeffs.num_row(), (unsigned int)2);
    ASSERT_EQ(testCoeffs.num_col(), (unsigned int)9);
    for (size_t i = 0; i < testCoeffs.num_row(); ++i) {
        for (size_t j = 0; j < testCoeffs.num_col(); ++j) {
            EXPECT_EQ(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1));
        }
    }
}

TEST(SquarePolynomialVectorTest, TestMakePolyVector) {
    OpalTestUtilities::SilenceTest silencer;

    std::vector<double> data(18, 0.);
    MMatrix<double> refCoeffs(2, 9, &data[0]);
    SquarePolynomialVector ref(2, refCoeffs);
    MVector<double> point(2, 0.);
    MVector<double> polyVector(9, -99);
    ref.MakePolyVector(point, polyVector);
    EXPECT_EQ(polyVector(1), 1.);
    for (size_t i = 1; i < 9; ++i)
        EXPECT_EQ(polyVector(i+1), 0.);
    point(1) = 3;
    point(2) = 2;
    //                    c   x0  x1 x0x1 x0x0 x0x0x1 x1x1 x1x1x0, x1x1x0x0
    double refVector[] = {1., 3., 2., 6., 9.,  18.,   4.,  12.,    36.};
    ref.MakePolyVector(point, polyVector);
    for (size_t i = 0; i < 9; ++i)
        EXPECT_EQ(polyVector(i+1), refVector[i]);
}

TEST(SquarePolynomialVectorTest, TestF) {
    OpalTestUtilities::SilenceTest silencer;

    std::vector<double> data(18);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;
    MMatrix<double> refCoeffs(2, 9, &data[0]);
    SquarePolynomialVector ref(2, refCoeffs);
    std::vector<double> point1(2, 0.);
    std::vector<double> value(2);
    ref.F(&point1[0], &value[0]);
    EXPECT_EQ(value[0], refCoeffs(1, 1));
    EXPECT_EQ(value[1], refCoeffs(2, 1));
    MVector<double> point2(2);
    point2(1) = 2.;
    point2(2) = 3.;
    MVector<double> polyVector(9, -99);
    ref.MakePolyVector(point2, polyVector);
    MVector<double> refValue = refCoeffs*polyVector;
    MVector<double> testValue(2, -99);
    ref.F(point2, testValue);
    EXPECT_EQ(testValue(1), refValue(1)); // sum (0, ..., 8)
    EXPECT_EQ(testValue(2), refValue(2)); // sum (9, ..., 17)
    std::vector<double> point3(2);
    point3[0] = 2.;
    point3[1] = 3.;
    ref.F(&point3[0], &value[0]);
    EXPECT_EQ(value[0], refValue(1)); // sum (0, ..., 8)
    EXPECT_EQ(value[1], refValue(2)); // sum (9, ..., 17)
}

TEST(SquarePolynomialVectorTest, TestDeriv) {
    std::vector<double> data(32);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;
    // Up to x^3
    MMatrix<double> inputCoeffs(2, 16, &data[0]);
    SquarePolynomialVector input(2, inputCoeffs);
    std::vector<int> derivPower({1, 2});
    SquarePolynomialVector test = input.Deriv(&derivPower[0]);
    std::vector<double> refCoeffs({
        28, 64, 156, 336, 132, 540, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        30, 68, 162, 348, 138, 558, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    });
    MMatrix<double> testCoeffs = test.GetCoefficientsAsMatrix();
    size_t index = 0;
    for (size_t i = 0; i < testCoeffs.num_row(); ++i) {
        for (size_t j = 0; j < testCoeffs.num_col(); ++j) {
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs[index], 1e-12);
            index++;
        }
    }
}
