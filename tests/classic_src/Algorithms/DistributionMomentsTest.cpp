//
// Class DistributionMoments
//   Computes the statistics of particle distributions.
//
// Copyright (c) 2021, Christof Metzger-Kraus
// All rights reserved
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

#include "gtest/gtest.h"

#include "Algorithms/DistributionMoments.h"
#include "Algorithms/OpalParticle.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"
#include "DistributionMomentsTestFixture.h"

TEST_F(DistributionMomentsTest, FullPositionTest) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedMeanPosition(5.75796e-8, 4.65360e-7, 1.56891e-3);
    EXPECT_NEAR(distributionMoments_m.getMeanPosition()(0), expectedMeanPosition(0), expectedMeanPosition(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getMeanPosition()(1), expectedMeanPosition(1), expectedMeanPosition(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getMeanPosition()(2), expectedMeanPosition(2), expectedMeanPosition(2) * 1e-4);

    Vector_t expectedStdPosition(1.85588e-4, 1.93126e-4, 2.68480e-4);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationPosition()(0), expectedStdPosition(0), expectedStdPosition(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationPosition()(1), expectedStdPosition(1), expectedStdPosition(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationPosition()(2), expectedStdPosition(2), expectedStdPosition(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, FullMomentumTest) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedMeanMomentum(1.47652e-6, 3.61031e-6, 2.44518e1);
    EXPECT_NEAR(distributionMoments_m.getMeanMomentum()(0), expectedMeanMomentum(0), expectedMeanMomentum(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getMeanMomentum()(1), expectedMeanMomentum(1), expectedMeanMomentum(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getMeanMomentum()(2), expectedMeanMomentum(2), expectedMeanMomentum(2) * 1e-4);

    Vector_t expectedStdMomentum(2.46592e-3, 2.49964e-3, 9.17394e-3);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationMomentum()(0), expectedStdMomentum(0), expectedStdMomentum(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationMomentum()(1), expectedStdMomentum(1), expectedStdMomentum(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getStandardDeviationMomentum()(2), expectedStdMomentum(2), expectedStdMomentum(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, FullNormalizeEmittanceTest) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedNormalizedEmittance(2.5302e-1, 2.55949e-1, 1.87916e-1);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance()(0) * 1e6,
                expectedNormalizedEmittance(0), expectedNormalizedEmittance(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance()(1) * 1e6,
                expectedNormalizedEmittance(1), expectedNormalizedEmittance(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance()(2) * 1e6,
                expectedNormalizedEmittance(2), expectedNormalizedEmittance(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, FullEnergyTest) {
    OpalTestUtilities::SilenceTest silencer;

    double expectedMeanGamma = 24.47225;
    EXPECT_NEAR(distributionMoments_m.getMeanGamma(),
                expectedMeanGamma, expectedMeanGamma * 1e-4);

    double expectedMeanKineticEnergy = 11.99427;
    EXPECT_NEAR(distributionMoments_m.getMeanKineticEnergy(),
                expectedMeanKineticEnergy, expectedMeanKineticEnergy * 1e-4);

    double expectedStdKineticEnergy = 4.68427e-3;
    EXPECT_NEAR(distributionMoments_m.getStdKineticEnergy(),
                expectedStdKineticEnergy, expectedStdKineticEnergy * 1e-4);

}

TEST_F(DistributionMomentsTest, FullDispersionTest) {
    OpalTestUtilities::SilenceTest silencer;

    double expectedDispersion_x = 1.41334e-6;
    EXPECT_NEAR(distributionMoments_m.getDx(),
                expectedDispersion_x, expectedDispersion_x * 1e-4);

    double expectedDispersion_px = 3.61358e-5;
    EXPECT_NEAR(distributionMoments_m.getDDx(),
                expectedDispersion_px, expectedDispersion_px * 1e-4);

    double expectedDispersion_y = 1.13848e-5;
    EXPECT_NEAR(distributionMoments_m.getDy(),
                expectedDispersion_y, expectedDispersion_y * 1e-4);

    double expectedDispersion_py = 8.83043e-5;
    EXPECT_NEAR(distributionMoments_m.getDDy(),
                expectedDispersion_py, expectedDispersion_py * 1e-4);
}

// To get the values for the percentiles in matlab extract the values from the
// fixture in the file DistributionMomentsTestFixture.cpp to e.g. a file
// distribution.txt. Then run

// A = load('distribution.txt');
// R = abs(A(:,1:2:5) - mean(A(:,1:2:5)));
// for d = 1:3
//   %percentile = prctile(R(:,d), 95.45);
//   [s I] = sort(R(:,d));
//   R = R(I,:);
//   N = round(size(A,1) * 0.9545);
//   percentile = (R(N) + R(N + 1)) / 2;
//   I = R(:,d) < percentile;
//   B = A(I, 2*d-1:2*d);
//   S = sqrt(1 - 1/N) * std(B);
//   M = mean(B);
//   RP = mean(B(:,1).*B(:,2)) - M(1)*M(2);
//   emittance = sqrt((S(1)*S(2))^2 - RP^2);
//   printf('%s: percentile = %g, emittance = %g\n', ...
//          ['X','Y','Z'](d), percentile, emittance);
// end

// Repeat the loop for the other percentiles.


TEST_F(DistributionMomentsTest, SixtyEightPercentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedPercentile(1.88262e-1, 1.96753e-1, 2.90999e-1);
    EXPECT_NEAR(distributionMoments_m.get68Percentile()(0) * 1e3,
                expectedPercentile(0), expectedPercentile(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get68Percentile()(1) * 1e3,
                expectedPercentile(1), expectedPercentile(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get68Percentile()(2) * 1e3,
                expectedPercentile(2), expectedPercentile(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NormalizedEmittanceAt68Percentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedEmittance(1.24194e-1, 1.23424e-1, 8.1483e-2);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance68Percentile()(0) * 1e6,
                expectedEmittance(0), expectedEmittance(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance68Percentile()(1) * 1e6,
                expectedEmittance(1), expectedEmittance(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance68Percentile()(2) * 1e6,
                expectedEmittance(2), expectedEmittance(2) * 1e-4);

}

TEST_F(DistributionMomentsTest, NinetyFivePercentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedPercentile(3.63691e-1, 3.77796e-1, 5.08294e-1);
    EXPECT_NEAR(distributionMoments_m.get95Percentile()(0) * 1e3,
                expectedPercentile(0), expectedPercentile(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get95Percentile()(1) * 1e3,
                expectedPercentile(1), expectedPercentile(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get95Percentile()(2) * 1e3,
                expectedPercentile(2), expectedPercentile(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NormalizedEmittanceAt95Percentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedEmittance(2.17939e-1, 2.20570e-1, 1.53668e-1);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance95Percentile()(0) * 1e6,
                expectedEmittance(0), expectedEmittance(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance95Percentile()(1) * 1e6,
                expectedEmittance(1), expectedEmittance(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance95Percentile()(2) * 1e6,
                expectedEmittance(2), expectedEmittance(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NinetyNinePercentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedPercentile(5.31749e-1, 5.56874e-1, 5.94257e-1);
    EXPECT_NEAR(distributionMoments_m.get99Percentile()(0) * 1e3,
                expectedPercentile(0), expectedPercentile(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get99Percentile()(1) * 1e3,
                expectedPercentile(1), expectedPercentile(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get99Percentile()(2) * 1e3,
                expectedPercentile(2), expectedPercentile(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NormalizedEmittanceAt99Percentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedEmittance(2.47584e-1, 2.50709e-1, 1.86971e-1);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99Percentile()(0) * 1e6,
                expectedEmittance(0), expectedEmittance(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99Percentile()(1) * 1e6,
                expectedEmittance(1), expectedEmittance(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99Percentile()(2) * 1e6,
                expectedEmittance(2), expectedEmittance(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NinetyNine_NinetyNinePercentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedPercentile(9.36452e-1, 9.83767e-1, 6.10324e-1);
    EXPECT_NEAR(distributionMoments_m.get99_99Percentile()(0) * 1e3,
                expectedPercentile(0), expectedPercentile(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get99_99Percentile()(1) * 1e3,
                expectedPercentile(1), expectedPercentile(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.get99_99Percentile()(2) * 1e3,
                expectedPercentile(2), expectedPercentile(2) * 1e-4);
}

TEST_F(DistributionMomentsTest, NormalizedEmittanceAt99_99Percentile) {
    OpalTestUtilities::SilenceTest silencer;

    Vector_t expectedEmittance(2.52308e-1, 2.55184e-1, 1.87881e-1);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99_99Percentile()(0) * 1e6,
                expectedEmittance(0), expectedEmittance(0) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99_99Percentile()(1) * 1e6,
                expectedEmittance(1), expectedEmittance(1) * 1e-4);
    EXPECT_NEAR(distributionMoments_m.getNormalizedEmittance99_99Percentile()(2) * 1e6,
                expectedEmittance(2), expectedEmittance(2) * 1e-4);
}