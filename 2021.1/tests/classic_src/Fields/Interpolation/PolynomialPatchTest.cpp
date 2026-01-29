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

#include "Fields/Interpolation/PolynomialPatch.h"

#include "opal_test_utilities/SilenceTest.h"

using namespace interpolation;

TEST(PolynomialPatchTest, TestPolynomialPatch) {
    OpalTestUtilities::SilenceTest silencer;

    ThreeDGrid grid(1, 2, 3, -1, -2, -3, 4, 3, 2);
   // we make a reference poly vector
    std::vector<double> data(54);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i/10.;
    MMatrix<double> refCoeffs(2, 27, &data[0]);
    SquarePolynomialVector ref(3, refCoeffs);
   // copy it into the grid
    std::vector<SquarePolynomialVector*> poly;
    for (int i = 0; i < grid.end().toInteger(); ++i)
        poly.push_back(new SquarePolynomialVector(ref));
    PolynomialPatch patch(grid.clone(), grid.clone(), poly);
    ThreeDGrid testGrid(1/4., 2/4., 3/4., -1, -2, -3, 4*4, 3*4, 2*4);
    for (Mesh::Iterator it = testGrid.begin(); it < testGrid.begin()+1; ++it) {
        std::vector<double> testValue(2);
        patch.function(&it.getPosition()[0], &testValue[0]);
        Mesh::Iterator nearest = grid.getNearest(&it.getPosition()[0]);
        std::vector<double> localPosition(3);
        for (size_t i = 0; i < 3; ++i) {
            localPosition[i] = nearest.getPosition()[i] - it.getPosition()[i];
        }
        std::vector<double> refValue(2);
        ref.F(&localPosition[0], &refValue[0]);
        for (size_t i = 0; i < 2; ++i)
            EXPECT_NEAR(testValue[i], refValue[i], 1e-6);
    }
}
