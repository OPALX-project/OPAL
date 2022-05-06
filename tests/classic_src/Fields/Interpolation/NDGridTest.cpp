/*
 *  Copyright (c) 2017, Chris Rogers
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
#include "Fields/Interpolation/NDGrid.h"
#include "Fields/Interpolation/Mesh.h"

#include "opal_test_utilities/SilenceTest.h"

namespace ndgridtest {

class NDGridTest : public ::testing::Test {
public:
    NDGridTest() : grid_m(nullptr) {
    }

    void SetUp( ) {
        std::vector< std::vector<double> > gridCoordinates(2);
        gridCoordinates[0] = std::vector<double>(2, 0.); // 0., 3.
        gridCoordinates[0][1] = 3.;
        gridCoordinates[1] = std::vector<double>(3, 1.); // 1., 5., 9.
        gridCoordinates[1][1] = 5.;
        gridCoordinates[1][2] = 9.;
        grid_m = new interpolation::NDGrid(gridCoordinates);
        gridCoordinates[1][2] = 10.; // force it to not be regular
        grid2_m = new interpolation::NDGrid(gridCoordinates);
    }

    void TearDown( ) {
        delete grid_m;
        grid_m = nullptr;
        delete grid2_m;
        grid2_m = nullptr;
    }

    ~NDGridTest() {
    }

    interpolation::NDGrid* grid_m;
    interpolation::NDGrid* grid2_m;
    OpalTestUtilities::SilenceTest silencer;

private:

};

TEST_F(NDGridTest, DefaultConstructorTest) {
    interpolation::NDGrid grid;
    EXPECT_EQ(grid.begin(), grid.end());
}

// also tests size
TEST_F(NDGridTest, Constructor1Test) {
    int size[] = {5, 6, 7, 8};
    double spacing[] = {1., 2., 3., 4.};
    double min[] = {-1., -2., -3., -4.};
    interpolation::NDGrid grid(4, size, spacing, min);

    ASSERT_EQ(grid.getPositionDimension(), 4);
    for (int i = 0; i < 4; ++i) {
        ASSERT_EQ(grid.size(i), size[i]);
        EXPECT_NEAR(grid.coord(1, i), min[i], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid.coord(size[i], i), min[i]+spacing[i]*(size[i]-1), 1e-12) << "Failed for i " << i;
    }
}

TEST_F(NDGridTest, Constructor2Test) {
    std::vector<int> size(2);
    size[0] = 2;
    size[1] = 3;
    std::vector<const double*> gridCoordinates(2, nullptr);
    const double vec1[] = {0., 2.};
    const double vec2[] = {1., 5., 9.};
    gridCoordinates[0] = vec1;
    gridCoordinates[1] = vec2;
    interpolation::NDGrid grid1(size, gridCoordinates);
    ASSERT_EQ(grid1.getPositionDimension(), 2);
    for (int i = 0; i < 2; ++i) {
        ASSERT_EQ(grid1.size(i), size[i]);
        EXPECT_NEAR(grid1.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid1.coord(i+2, i), gridCoordinates[i][size[i]-1], 1e-12) << "Failed for i " << i;
    }
}

TEST_F(NDGridTest, Constructor3Test) {
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;
    interpolation::NDGrid grid1(gridCoordinates);
    ASSERT_EQ(grid1.getPositionDimension(), 2);
    for (int i = 0; i < 2; ++i) {
        int size = gridCoordinates[i].size();
        ASSERT_EQ(grid1.size(i), size);
        EXPECT_NEAR(grid1.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid1.coord(i+2, i), gridCoordinates[i][size-1], 1e-12) << "Failed for i " << i;
    }
}

TEST_F(NDGridTest, CpyConstructorTest) {
    interpolation::NDGrid grid(*grid_m);
    EXPECT_EQ(grid.getConstantSpacing(), grid_m->getConstantSpacing());
    ASSERT_EQ(grid.getPositionDimension(), grid_m->getPositionDimension());
    for (int i = 0; i < grid.getPositionDimension(); ++i) {
        ASSERT_EQ(grid.coordVector(i).size(), grid_m->coordVector(i).size());
        for (size_t j = 0; j < grid.coordVector(i).size(); ++j) {
            EXPECT_NEAR(grid.coordVector(i)[j],
                        grid_m->coordVector(i)[j], 1e-15);
        }
    }
}

TEST_F(NDGridTest, CoordTest) {
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;
    interpolation::NDGrid grid_var(gridCoordinates);
    const interpolation::NDGrid grid_const(gridCoordinates);
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(grid_var.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid_const.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
    }
}

TEST_F(NDGridTest, CoordVectorTest) {  // and newCoordArray
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;
    interpolation::NDGrid grid(gridCoordinates);
    for (int i = 0; i < 2; ++i) {
        std::vector<double> coords_v = grid.coordVector(i);
        double* coords_a = grid.newCoordArray(i);
        ASSERT_EQ(coords_v.size(), gridCoordinates[i].size());
        for (size_t j = 0; j < gridCoordinates[i].size(); j++) {
            EXPECT_NEAR(coords_v[j], gridCoordinates[i][j], 1e-12);
            EXPECT_NEAR(coords_a[j], gridCoordinates[i][j], 1e-12);
        }
        delete[] coords_a;
    }
}

TEST_F(NDGridTest, CoordLowerBoundTest) {
    // first dimension ... 0., 3.;
    int index = -1;
    grid_m->coordLowerBound(-1., 0, index);
    EXPECT_EQ(index, -1);
    index = -2;
    grid_m->coordLowerBound(0.0001, 0, index);
    EXPECT_EQ(index, 0);
    index = -2;
    grid_m->coordLowerBound(2.999, 0, index);
    EXPECT_EQ(index, 0);
    index = -2;
    grid_m->coordLowerBound(3.001, 0, index);
    EXPECT_EQ(index, 1);
    index = -2;
    grid_m->coordLowerBound(1000.0001, 0, index);
    EXPECT_EQ(index, 1);

    // second dimension ...  1., 5., 9.
    index = -2;
    grid_m->coordLowerBound(0.999, 1, index);
    EXPECT_EQ(index, -1);
    index = -2;
    grid_m->coordLowerBound(4.999, 1, index);
    EXPECT_EQ(index, 0);
    index = -2;
    grid_m->coordLowerBound(8.999, 1, index);
    EXPECT_EQ(index, 1);
    index = -2;
    grid_m->coordLowerBound(9.0001, 1, index);
    EXPECT_EQ(index, 2);
    index = -2;
    grid_m->coordLowerBound(1000.0001, 1, index);
    EXPECT_EQ(index, 2);
}

TEST_F(NDGridTest, LowerBoundTest) {
    // first dimension ... 0., 3.;
    // second dimension ...  1., 5., 9.
    std::vector<int> index1(2, -2);
    std::vector<double> pos1(2, -1.);
    grid_m->lowerBound(pos1, index1);
    EXPECT_EQ(index1[0], -1);
    EXPECT_EQ(index1[1], -1);

    std::vector<int> index2(2, -2);
    std::vector<double> pos2(2, 4.);
    grid_m->lowerBound(pos2, index2);
    EXPECT_EQ(index2[0], 1);
    EXPECT_EQ(index2[1], 0);

    std::vector<int> index3(2, -2);
    std::vector<double> pos3(2, 100.);
    grid_m->lowerBound(pos3, index3);
    EXPECT_EQ(index3[0], 1);
    EXPECT_EQ(index3[1], 2);
}

TEST_F(NDGridTest, MinMaxTest) {
    EXPECT_EQ(grid_m->min(0), 0.);
    EXPECT_EQ(grid_m->min(1), 1.);
    EXPECT_EQ(grid_m->max(0), 3.);
    EXPECT_EQ(grid_m->max(1), 9.);
}

TEST_F(NDGridTest, SetCoordTest) {
    double xNew[] = {5., 10., 12., 15.};
    grid_m->setCoord(0, 4, xNew);
    std::vector<double> xTest = grid_m->coordVector(0);
    EXPECT_EQ(xTest.size(), (unsigned int)4);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(xTest[i], xNew[i]);
    }
}

TEST_F(NDGridTest, BeginEndTest) {
    ASSERT_EQ(grid_m->begin().getState().size(), (unsigned int)2);
    EXPECT_EQ(grid_m->begin().getState()[0], 1);
    EXPECT_EQ(grid_m->begin().getState()[1], 1);
    ASSERT_EQ(grid_m->end().getState().size(), (unsigned int)2);
    EXPECT_EQ(grid_m->end().getState()[0], 3); // one past the last (2, 3)
    EXPECT_EQ(grid_m->end().getState()[1], 1);
    EXPECT_EQ(grid_m->begin().toInteger(), 0);
    EXPECT_EQ(grid_m->end().toInteger(), 6);
    
}

TEST_F(NDGridTest, GetPositionTest) {
    std::vector<double> position(3, -1);
    interpolation::Mesh::Iterator it = grid_m->begin();
    grid_m->getPosition(it, &position[0]);
    EXPECT_EQ(grid_m->getPositionDimension(), 2);
    EXPECT_EQ(position[0], 0.);
    EXPECT_EQ(position[1], 1.);
    it[0] = 2;
    it[1] = 3;
    grid_m->getPosition(it, &position[0]);
    EXPECT_EQ(position[0], 3.);
    EXPECT_EQ(position[1], 9.);
}

TEST_F(NDGridTest, GetSetConstantSpacingTest) {
    EXPECT_TRUE(grid_m->getConstantSpacing());
    grid_m->setConstantSpacing(false);
    EXPECT_FALSE(grid_m->getConstantSpacing());
    grid_m->setConstantSpacing(true);
    EXPECT_TRUE(grid_m->getConstantSpacing());
    grid_m->setConstantSpacing(false);
    grid_m->setConstantSpacing();
    EXPECT_TRUE(grid_m->getConstantSpacing());

    EXPECT_FALSE(grid2_m->getConstantSpacing());
    grid2_m->setConstantSpacing(true);
    EXPECT_TRUE(grid2_m->getConstantSpacing());
    grid2_m->setConstantSpacing(false);
    EXPECT_FALSE(grid2_m->getConstantSpacing());
    grid2_m->setConstantSpacing(true);
    grid2_m->setConstantSpacing(0.1);
    EXPECT_FALSE(grid2_m->getConstantSpacing());
    grid2_m->setConstantSpacing(2.01);
    EXPECT_TRUE(grid2_m->getConstantSpacing());
}

TEST_F(NDGridTest, ToIntegerTest) {
    interpolation::Mesh::Iterator it = grid_m->begin();
    EXPECT_EQ(grid_m->toInteger(it), 0);
    it[0] = 2;
    it[1] = 3;
    EXPECT_EQ(grid_m->toInteger(it), 5);
}

TEST_F(NDGridTest, GetNearestTest) {

    // first dimension ... 0., 3.;
    // second dimension ...  1., 5., 9.

    std::vector<double> pos(2, 0.);
    interpolation::Mesh::Iterator it;
    pos[0] = 1.49;
    pos[1] = 2.99;
    it = grid_m->getNearest(&pos[0]);
    EXPECT_EQ(it[0], 1);
    EXPECT_EQ(it[1], 1);
    pos[0] = 0.49;
    it = grid_m->getNearest(&pos[0]);
    EXPECT_EQ(it[0], 1);
    EXPECT_EQ(it[1], 1);
    pos[1] = 9.49;
    it = grid_m->getNearest(&pos[0]);
    EXPECT_EQ(it[0], 1);
    EXPECT_EQ(it[1], 3);

    pos[0] = 3.49;
    it = grid_m->getNearest(&pos[0]);
    EXPECT_EQ(it[0], 2);
    EXPECT_EQ(it[1], 3);
}

TEST_F(NDGridTest, DualTest) {
    using interpolation::NDGrid;
    NDGrid* new_grid = dynamic_cast<NDGrid*>(grid2_m->dual());
    std::vector<std::vector<double> > ref = {
      {1.5},
      {3.0, 7.5}
    };
    ASSERT_EQ(new_grid->getPositionDimension(), int(ref.size()));
    for (size_t i = 0; i < ref.size(); ++i) {
        ASSERT_EQ(new_grid->size(i), int(ref[i].size()));
        for (size_t j = 0; j < ref[i].size(); ++j) {
            EXPECT_NEAR(new_grid->coordVector(i)[j], ref[i][j], 1e-9) << i << " " << j;
        }
    }
}

TEST_F(NDGridTest, IsOutOfBoundsBeginTest) {
    // nb isOutOfBounds is defined in Mesh.hh (but tested here for convenience
    interpolation::Mesh::Iterator it = grid_m->begin();
    EXPECT_FALSE(it.isOutOfBounds());
    for (size_t i = 0; i < it.getState().size(); ++i) {
        it[i] -= 1;
        EXPECT_TRUE(it.isOutOfBounds());
        it[i] += 1;
    }
}

TEST_F(NDGridTest, IsOutOfBoundsEndTest) {
    interpolation::Mesh::Iterator it = grid_m->end();
    it--;
    EXPECT_FALSE(it.isOutOfBounds());
    for (size_t i = 0; i < it.getState().size(); ++i) {
        it[i] += 1;
        EXPECT_TRUE(it.isOutOfBounds());
        it[i] -= 1;
    }
}

} // namespace ndgridtest
