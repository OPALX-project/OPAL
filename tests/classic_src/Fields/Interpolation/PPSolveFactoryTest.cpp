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

#include <sstream>

#include "gtest/gtest.h"

//  #define PolynomialPatchTest_MakePlots
#ifdef PolynomialPatchTest_MakePlots
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#endif

#include "Fields/Interpolation/SquarePolynomialVector.h"
#include "Fields/Interpolation/PolynomialPatch.h"
#include "Fields/Interpolation/PolynomialCoefficient.h"
#include "Fields/Interpolation/PPSolveFactory.h"

#include "opal_test_utilities/SilenceTest.h"

using namespace interpolation;

namespace {
void test_points(int dim, int lower, int upper, std::vector< std::vector<int> > pts) {
      int upper_size = 1;
      int lower_size = 1;
      for (int i = 0; i < dim; ++i)
          upper_size *= upper+1;
      for (int i = 0; i < dim; ++i)
          lower_size *= lower+1;
      if (lower < 0)
          lower_size = 0;
      if (upper < 0)
          upper_size = 0;
      // size should be difference in area of the squares
      EXPECT_EQ(pts.size(), (unsigned int)(upper_size - lower_size));
      for (size_t i = 0; i < pts.size(); ++i) {
          // each pts element should have length dim
          EXPECT_EQ(pts[i].size(), (size_t)dim);
          // each pts element should have indices with lower < size <= upper
          bool in_bounds = true;
          for (int j = 0; j < dim; ++j) {
              in_bounds = in_bounds || (pts[i][j] > lower && pts[i][j] <= upper);
          }
          EXPECT_TRUE(in_bounds);
          // each pts element should be unique
          for (size_t j = 0; j < pts.size(); ++j) {
              if (j == i)
                  continue;
              bool equal = true;
              for (int k = 0; k < dim; ++k)
                  equal &= pts[i][k] == pts[j][k];
              EXPECT_FALSE(equal);
          }
      }
}
}

TEST(PPSolveFactoryTest, TestNearbyPointsSquares) {
    OpalTestUtilities::SilenceTest silencer;

    for (int upper = 0; upper < 5; ++upper) {
        for (int lower = 0; lower < upper; ++lower) {
            for (int dim = 1; dim < 5; ++dim) {
                std::vector< std::vector<int> > pts =
                    PPSolveFactory::getNearbyPointsSquares(dim, lower, upper);
                test_points(dim, lower, upper, pts);
            }
        }
    }
}

class PPSolveFactoryTestFixture : public ::testing::Test {
  protected:
    std::vector<std::vector<double> > values;
    SquarePolynomialVector ref;
    ThreeDGrid* grid;
    MMatrix<double> refCoeffs;
    int np;

    PPSolveFactoryTestFixture() {
       // we make a reference poly vector
        std::vector<double> data(2*27); // 2^3
        for (size_t i = 0; i < data.size(); ++i)
            data[i] = i/10.;
        refCoeffs = MMatrix<double>(2, 27, &data[0]);
        ref = SquarePolynomialVector(3, refCoeffs);
        // we get some data
        np = 6;
        // choose the grid so that a midpoint falls at 0/0/0; we want to compare it
        // later on...
        grid = new ThreeDGrid(1., 2., 3., -2.5, -5.0, -7.5, np, np, np);
        for (Mesh::Iterator it = grid->begin(); it < grid->end(); ++it) {
            std::vector<double> value(2);
            ref.F(&it.getPosition()[0], &value[0]);
            values.push_back(value);
        }
    }
};

// check linear fit exactly reproduces data on grid points
TEST_F(PPSolveFactoryTestFixture, TestSolvePolynomialLinear) {
    OpalTestUtilities::SilenceTest silencer;

    PPSolveFactory fac1(grid->clone(),
                       values,
                       1,
                       1);
    PolynomialPatch* patch1 = fac1.solve();
    for (Mesh::Iterator it = grid->begin(); it < grid->end(); ++it) {
        std::vector<double> value(2);
        patch1->function(&it.getPosition()[0], &value[0]);
        for (size_t i = 0; i < 2; ++i)
            EXPECT_NEAR(values[it.toInteger()][i], value[i], 1e-6);
    }
    delete patch1;
}

// check quadratic fit exactly reproduces data on and off grid points (data
// comes from a quadratic polynomial as source)
TEST_F(PPSolveFactoryTestFixture, TestSolvePolynomialQuadratic) {
    OpalTestUtilities::SilenceTest silencer;

    PPSolveFactory fac2(grid->clone(),
                       values,
                       2,
                       2);
    PolynomialPatch* patch2 = fac2.solve();
    // first check that the polynomial at 0, 0, 0 is the same as ref
    std::vector<double> zero(3, 0.);
    SquarePolynomialVector* test = patch2->getPolynomialVector(&zero[0]);
    EXPECT_EQ(patch2->getValueDimension(), (unsigned int)2);
    MMatrix<double> testCoeffs = test->GetCoefficientsAsMatrix();
    std::cout << "Ref" << std::endl;
    std::cout << refCoeffs << std::endl;
    std::cout << "Test" << std::endl;
    std::cout << testCoeffs << std::endl;
    ASSERT_EQ(testCoeffs.num_row(), refCoeffs.num_row());
    ASSERT_EQ(testCoeffs.num_col(), refCoeffs.num_col());
    for (size_t i = 0; i < testCoeffs.num_row(); ++i)
        for (size_t j = 0; j < testCoeffs.num_row(); ++j) {
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6)
                                               << "col " << i << " row " << j;
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6)
                                               << "col " << i << " row " << j;
        }

    // now check that the values in each polynomial are correct
    ThreeDGrid testGrid(1./4., 2./4., 3./4., -1., -2., -3., np*4-1, np*4-1, np*4-1);
    for (Mesh::Iterator it = testGrid.begin(); it < testGrid.end(); ++it) {
        std::vector<double> refValue(2);
        std::vector<double> testValue(2);
        ref.F(&it.getPosition()[0], &refValue[0]);
        patch2->function(&it.getPosition()[0], &testValue[0]);
        for (size_t i = 0; i < 2; ++i)
            EXPECT_NEAR(refValue[i], testValue[i], 1e-6) << std::endl << it;
    }
}

// check smoothed quadratic fit exactly reproduces data on and off grid points
// except near to the boundary
TEST_F(PPSolveFactoryTestFixture, DISABLED_TestSolvePolynomialQuadraticSmoothed) {
    OpalTestUtilities::SilenceTest silencer;

    PPSolveFactory fac2(grid->clone(),
                       values,
                       1,
                       2);
    PolynomialPatch* patch2 = fac2.solve();
    // first check that the polynomial at 0, 0, 0 is the same as ref
    std::vector<double> zero(3, 0.);
    SquarePolynomialVector* test = patch2->getPolynomialVector(&zero[0]);
    MMatrix<double> testCoeffs = test->GetCoefficientsAsMatrix();
    std::cout << "Ref" << std::endl;
    std::cout << refCoeffs << std::endl;
    std::cout << "Test" << std::endl;
    std::cout << testCoeffs << std::endl;
    ASSERT_EQ(testCoeffs.num_row(), refCoeffs.num_row());
    ASSERT_EQ(testCoeffs.num_col(), refCoeffs.num_col());
    for (size_t i = 0; i < testCoeffs.num_row(); ++i)
        for (size_t j = 0; j < testCoeffs.num_row(); ++j) {
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6)
                                               << "col " << i << " row " << j;
            EXPECT_NEAR(testCoeffs(i+1, j+1), refCoeffs(i+1, j+1), 1e-6)
                                               << "col " << i << " row " << j;
        }

    // now check that the values in each polynomial are correct
    ThreeDGrid testGrid(1./4., 2./4., 3./4., -1., -2., -3., np*4-1, np*4-1, np*4-1);
    for (Mesh::Iterator it = testGrid.begin(); it < testGrid.end(); ++it) {
        std::vector<double> refValue(2);
        std::vector<double> testValue(2);
        ref.F(&it.getPosition()[0], &refValue[0]);
        patch2->function(&it.getPosition()[0], &testValue[0]);
        for (size_t i = 0; i < 2; ++i)
            EXPECT_NEAR(refValue[i], testValue[i], 1e-6) << std::endl << it;
    }
}

namespace {
std::vector<double> get_value(std::vector<double> x, int np) {
    static const double twopi = asin(1.)*4;
    std::vector<double> a_value(3);
    a_value[0] = sin(twopi*x[0]/np)*sin(twopi*x[1]/np)*sin(twopi*x[2]/np);
    a_value[1] = cos(twopi*x[0]/np)*cos(twopi*x[1]/np)*cos(twopi*x[2]/np);
    a_value[2] = 1.;
    return a_value;
}

#ifdef PolynomialPatchTest_MakePlots
void plot(int n_points, std::vector<double> start, std::vector<double> end, PolynomialPatch* patch, int n_grid_points, std::string title) {
    double delta_sq = 0;
    for (size_t i = 0; i < start.size(); ++i)
        delta_sq += (start[i]-end[i])*(start[i]-end[i]);
    std::vector<TGraph*> graphs(3*patch->getValueDimension(), NULL);
    for (size_t j = 0; j < graphs.size(); ++j) {
        graphs[j] = new TGraph(n_points);
    }
    std::vector<double> errors(3, 0.);
    for (int i = 0; i < n_points; ++i) {
        std::vector<double> point = start;
        for (size_t j = 0; j < point.size(); ++j) {
            point[j] += (end[j]-start[j])*i/n_points;
        }
        std::vector<double> fit_value(patch->getValueDimension(), 0.);
        patch->function(&point[0], &fit_value[0]);
        for (size_t j = 0; j < patch->getValueDimension(); ++j) {
            graphs[j]->SetPoint(i, point[0], fit_value[j]);
        }
        std::vector<double> true_value = get_value(point, n_grid_points);
        for (size_t j = 0; j < patch->getValueDimension(); ++j) {
            graphs[j+patch->getValueDimension()]->SetPoint(i, point[0], true_value[j]);
            graphs[j+2*patch->getValueDimension()]->SetPoint(i, point[0], fabs(fit_value[j]-true_value[j]));
            errors[j] += fabs(fit_value[j]-true_value[j]);
        }
    }
    graphs[0]->SetTitle("Fit Bx");
    graphs[1]->SetTitle("Fit By");
    graphs[2]->SetTitle("Fit Bz");
    graphs[3]->SetTitle("True Bx");
    graphs[4]->SetTitle("True By");
    graphs[5]->SetTitle("True Bz");

    TCanvas* c1 = new TCanvas("canvas", "canvas", 1024, 768);
    c1->Draw();
    TH2D* h1 = new TH2D("h1", (title+";pos [AU]; Fit [AU]").c_str(), 10000, start[0]-(end[0]-start[0])/10., end[0]+(end[0]-start[0])/10., 1000, -5., 5.);
    h1->SetStats(false);
    h1->Draw();
    for (size_t j = 0; j < 2*patch->getValueDimension(); ++j) {
        graphs[j]->SetFillColor(0);
        graphs[j]->SetLineColor(j+1);
        graphs[j]->SetMarkerColor(j+1);
        graphs[j]->Draw("l");
    }
    c1->BuildLegend();
    static int test_index = 0;
    std::stringstream test_name;
    test_name << "test_" << test_index;
    c1->SetGridx();
    c1->Print((test_name.str()+".png").c_str());
    c1->Print((test_name.str()+".root").c_str());

    TCanvas* c2 = new TCanvas("deltas canvas", "deltas canvas", 1024, 768);
    c2->Draw();
    TH2D* h2 = new TH2D("h2", ("Residuals "+title+";pos [AU]; Fit [AU]").c_str(), 10000, start[0]-(end[0]-start[0])/10., end[0]+(end[0]-start[0])/10., 1000, 1e-6, 5.);
    h2->SetStats(false);
    h2->Draw();

    std::vector<std::stringstream> strings(errors.size());
    for (size_t i = 0; i < errors.size(); ++i) {
        strings[i] << " <#epsilon> " << errors[i]/n_points;
    }

    graphs[6]->SetTitle(("Delta Bx "+strings[0].str()).c_str());
    graphs[7]->SetTitle(("Delta By "+strings[1].str()).c_str());
    graphs[8]->SetTitle(("Delta Bz "+strings[2].str()).c_str());
    for (size_t j = 2*patch->getValueDimension(); j < 3*patch->getValueDimension(); ++j) {
        graphs[j]->SetFillColor(0);
        graphs[j]->SetLineColor(j+1);
        graphs[j]->SetMarkerColor(j+1);
        graphs[j]->Draw("l");
    }
    c2->BuildLegend();
    c2->SetGridx();
    c2->SetLogy();
    c2->Print(("residuals_"+test_name.str()+".png").c_str());
    c2->Print(("residuals_"+test_name.str()+".root").c_str());
    test_index++;
}
#else // PolynomialPatchTest_MakePlots
    void plot(int /*n_points*/, std::vector<double> /*start*/, std::vector<double> /*end*/, PolynomialPatch* /*patch*/, int /*n_grid_points*/, std::string /*title*/) {
}
#endif // PolynomialPatchTest_MakePlots
}

TEST(PPSolveFactoryTest, TestThreeDSolveSinCos) {
    OpalTestUtilities::SilenceTest silencer;

    int np = 11;
    ThreeDGrid grid(1., 1., 1., 1., 2., 3., np, np, np);
    // test grid starts a few points in to avoid boundary effects
    ThreeDGrid fineGrid(1./4., 1./4., 1./4., 3., 4., 5., (np-2)*4, (np-2)*4, (np-2)*4);
    std::vector<double> start = grid.begin().getPosition();
    std::vector<double> end = (grid.end()-1).getPosition();

    std::vector<std::vector<double> > values;
    for (Mesh::Iterator it = grid.begin(); it < grid.end(); ++it) {
        values.push_back(get_value(it.getPosition(), np));
    }
    std::vector<PolynomialPatch*> ppVec;
    for (int smooth_order = 1; smooth_order < 4; ++smooth_order) {
        for (int pp_order = smooth_order; pp_order <= smooth_order; ++pp_order) {
            std::cout << "Building pp of order " << pp_order << " smooth " << smooth_order << std::endl;
            PPSolveFactory fac(grid.clone(),
                               values,
                               pp_order,
                               smooth_order);
            ppVec.push_back(fac.solve());
            std::stringstream ss;
            ss << "smooth: " << smooth_order << " poly: " << pp_order;
            plot(100, start, end, ppVec.back(), np, ss.str());
            EXPECT_EQ(ppVec.back()->getValueDimension(), (unsigned int)3);
        }
    }
    std::cout << "Testing" << std::endl;
    // check we get exactly right answer on mesh points
    for (Mesh::Iterator it = grid.begin(); it < grid.end(); ++it) {
        for (size_t ppi = 0; ppi < ppVec.size(); ppi++) {
            std::vector<double> testValue(3, 0.);
            ppVec[ppi]->function(&it.getPosition()[0], &testValue[0]);
            for (size_t i = 0; i < 3; ++i) {
                EXPECT_NEAR(get_value(it.getPosition(), np)[i], testValue[i], 1e-6);
            }
        }
    }
    // check we get nearly right answer between mesh points
    std::vector<std::vector<double> > sumError(ppVec.size(), std::vector<double>(3, 0.));
    for (Mesh::Iterator it = fineGrid.begin(); it < fineGrid.end(); ++it) {
        std::vector<double> prevTestValue(3, 0.);
        std::vector<double> refValue = get_value(it.getPosition(), np);
        ppVec[0]->function(&it.getPosition()[0], &prevTestValue[0]);
        for (size_t ppi = 0; ppi < ppVec.size(); ppi++) {
            std::vector<double> testValue(3, 0.);
            ppVec[ppi]->function(&it.getPosition()[0], &testValue[0]);
            for (size_t i = 0; i < 3; ++i) {
                sumError[ppi][i] += fabs(testValue[i]-refValue[i]);
            }
            prevTestValue = testValue;
        }
    }

    for (size_t i = 0; i < sumError.size(); ++i) {
        std::cout << i << ": ";
        for (size_t j = 0; j < sumError[i].size(); ++j) {
            std::cout << sumError[i][j] << " ";
            // higher order interpolations should fit better than lower order
            // except a little bit of floating point error can come in for e.g.
            // Bz where we are fitting a straight line (pathological case)...
            if (i > 0) {
                EXPECT_LT(sumError[i][j], sumError[i-1][j]+1e-6);
            }
        }
        std::cout << std::endl;
    }
}