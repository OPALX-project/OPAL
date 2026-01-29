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
#include "Fields/Interpolation/ThreeDGrid.h"

#include "opal_test_utilities/SilenceTest.h"

using interpolation::ThreeDGrid;

namespace threegridtest {

TEST(ThreeDGridTest, LowerBoundTest) {
    OpalTestUtilities::SilenceTest silencer;

    std::vector<double> xVar(7), yVar(8), zVar(9);
    for (size_t i = 0; i < xVar.size(); ++i) {
        xVar[i] = 4.+1.*i;
    }
    for (size_t i = 0; i < yVar.size(); ++i) {
        yVar[i] = -5.+2.*i;
    }
    for (size_t i = 0; i < zVar.size(); ++i) {
        zVar[i] = 6.+3.*i*i;
    }
    zVar[5] = (zVar[4]+zVar[6])*0.51; // force not constant spacing
    ThreeDGrid gridConst(1., 2., 3., // dx/dy/dz
                    4., 5., 6., // min x/y/z
                    7, 8, 9); // n x/y/z coordinates in the grid
    ThreeDGrid gridVar(xVar, yVar, zVar);
    gridVar.setConstantSpacing(false);

    ThreeDGrid* gridArray[] = {&gridConst, &gridVar};
    for (size_t i = 0; i < 2; ++i) {
        ThreeDGrid* grid = gridArray[i];
        int index = -99;
        for (int j = 0; j < grid->xSize(); ++j) {
            grid->xLowerBound(grid->x(j+1)-1e-9, index);
            EXPECT_EQ(index, int(j)-1) << "grid" << i << " " << j;
            grid->xLowerBound(grid->x(j+1)+0e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
            grid->xLowerBound(grid->x(j+1)+1e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
        }

        index = -99;
        for (int j = 0; j < grid->ySize(); ++j) {
            grid->yLowerBound(grid->y(j+1)-1e-9, index);
            EXPECT_EQ(index, int(j)-1) << "grid" << i << " " << j;
            grid->yLowerBound(grid->y(j+1)+0e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
            grid->yLowerBound(grid->y(j+1)+1e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
        }

        index = -99;
        for (int j = 0; j < grid->zSize(); ++j) {
            grid->zLowerBound(grid->z(j+1)-1e-9, index);
            EXPECT_EQ(index, int(j)-1) << "grid" << i << " " << j;
            grid->zLowerBound(grid->z(j+1)+0e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
            grid->zLowerBound(grid->z(j+1)+1e-9, index);
            EXPECT_EQ(index, int(j)) << "grid" << i << " " << j;
        }
    }
}
}
