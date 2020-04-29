#include "gtest/gtest.h"

#include "Distribution/Distribution.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"

#include "gsl/gsl_statistics_double.h"

TEST(GaussTest, FullSigmaTest1) {
    OpalTestUtilities::SilenceTest silencer;

    const double expectedR11 = 1.978;
    const double expectedR22 = 0.7998;
    const double expectedR33 = 2.498;
    const double expectedR44 = 0.6212;
    const double expectedR55 = 1.537;
    const double expectedR66 = 0.9457;

    const double expectedR21 = -0.40993;
    const double expectedR43 = 0.77208;
    const double expectedR65 = 0.12051;
    const double expectedR51 = 0.14935;
    const double expectedR52 = 0.59095;
    const double expectedR61 = 0.72795;
    const double expectedR62 = -0.3550;

    std::vector<double> expectedR({expectedR21, 0, 0,           expectedR51, expectedR61, \
                /*                           */ 0, 0,           expectedR52, expectedR62, \
                /*                              */ expectedR43, 0,           0, \
                /*                                            */0,           0,                                         \
                /*                                                         */expectedR65});

    Distribution dist;

    Attributes::setUpperCaseString(dist.itsAttr[Attrib::Distribution::TYPE], "GAUSS");
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAX], expectedR11 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPX], expectedR22);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAY], expectedR33 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPY], expectedR44);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAZ], expectedR55 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPZ], expectedR66);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFX], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFY], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFLONG], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPX], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPY], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPZ], 5.0);
    Attributes::setRealArray(dist.itsAttr[Attrib::Distribution::R], expectedR);
    Attributes::setBool(dist.itsAttr[Attrib::Distribution::EMITTED], false);
    Attributes::setBool(dist.itsAttr[Attrib::Distribution::WRITETOFILE], true);

    dist.setDistType();
    dist.checkIfEmitted();
    size_t numParticles = 1000000;
    dist.totalNumberParticles_m = numParticles;
    dist.create(numParticles, Physics::m_p);

    double R11 = sqrt(gsl_stats_variance(&(dist.xDist_m[0]), 1, dist.xDist_m.size())) * 1e3;
    double R22 = sqrt(gsl_stats_variance(&(dist.pxDist_m[0]), 1, dist.pxDist_m.size()));
    double R33 = sqrt(gsl_stats_variance(&(dist.yDist_m[0]), 1, dist.yDist_m.size())) * 1e3;
    double R44 = sqrt(gsl_stats_variance(&(dist.pyDist_m[0]), 1, dist.pyDist_m.size()));
    double R55 = sqrt(gsl_stats_variance(&(dist.tOrZDist_m[0]), 1, dist.tOrZDist_m.size())) * 1e3;
    double R66 = sqrt(gsl_stats_variance(&(dist.pzDist_m[0]), 1, dist.pzDist_m.size()));

    double R21 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.pxDist_m[0]), 1, dist.xDist_m.size()) * 1e3 /
                  (expectedR11 * expectedR22));
    double R43 = (gsl_stats_covariance(&(dist.yDist_m[0]), 1, &(dist.pyDist_m[0]), 1, dist.yDist_m.size()) * 1e3 /
                  (expectedR33 * expectedR44));
    double R51 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.tOrZDist_m[0]), 1, dist.xDist_m.size()) * 1e6 /
                  (expectedR11 * expectedR55));
    double R52 = (gsl_stats_covariance(&(dist.pxDist_m[0]), 1, &(dist.tOrZDist_m[0]), 1, dist.pxDist_m.size()) * 1e3 /
                  (expectedR22 * expectedR55));
    double R61 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.pzDist_m[0]), 1, dist.xDist_m.size()) * 1e3 /
                  (expectedR11 * expectedR66));
    double R62 = (gsl_stats_covariance(&(dist.pxDist_m[0]), 1, &(dist.pzDist_m[0]), 1, dist.pxDist_m.size()) /
                  (expectedR22 * expectedR66));

    EXPECT_NEAR(R11 / expectedR11, 1.0, 0.002);
    EXPECT_NEAR(R22 / expectedR22, 1.0, 0.002);
    EXPECT_NEAR(R33 / expectedR33, 1.0, 0.002);
    EXPECT_NEAR(R44 / expectedR44, 1.0, 0.002);
    EXPECT_NEAR(R55 / expectedR55, 1.0, 0.002);
    EXPECT_NEAR(R66 / expectedR66, 1.0, 0.002);

    EXPECT_NEAR(R21 / expectedR21, 1.0, 0.02);
    EXPECT_NEAR(R43 / expectedR43, 1.0, 0.02);
    EXPECT_NEAR(R51 / expectedR51, 1.0, 0.05);
    EXPECT_NEAR(R52 / expectedR52, 1.0, 0.02);
    EXPECT_NEAR(R61 / expectedR61, 1.0, 0.02);
    EXPECT_NEAR(R62 / expectedR62, 1.0, 0.02);
}

TEST(GaussTest, FullSigmaTest2) {
    OpalTestUtilities::SilenceTest silencer;

    const double expectedR11 = 1.978;
    const double expectedR22 = 0.7998;
    const double expectedR33 = 2.498;
    const double expectedR44 = 0.6212;
    const double expectedR55 =  1.537;
    const double expectedR66 = 0.9457;

    const double expectedR21 = -0.40993;
    const double expectedR43 = 0.77208;
    const double expectedR65 = 0.12051;
    const double expectedR51 = 0.14935;
    const double expectedR52 = 0.59095;
    const double expectedR61 = 0.72795;
    const double expectedR62 = -0.3550;

    Distribution dist;

    Attributes::setUpperCaseString(dist.itsAttr[Attrib::Distribution::TYPE], "GAUSS");
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAX], expectedR11 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPX], expectedR22);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAY], expectedR33 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPY], expectedR44);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAZ], expectedR55 * 1e-3);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::SIGMAPZ], expectedR66);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CORRX], expectedR21);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CORRY], expectedR43);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CORRZ], expectedR65);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::R51], expectedR51);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::R61], expectedR61);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::R52], expectedR52);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::R62], expectedR62);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFX], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFY], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFLONG], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPX], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPY], 5.0);
    Attributes::setReal(dist.itsAttr[Attrib::Distribution::CUTOFFPZ], 5.0);
    Attributes::setBool(dist.itsAttr[Attrib::Distribution::EMITTED], false);
    Attributes::setBool(dist.itsAttr[Attrib::Distribution::WRITETOFILE], true);

    dist.setDistType();
    dist.checkIfEmitted();

    size_t numParticles = 1000000;
    dist.totalNumberParticles_m = numParticles;
    dist.create(numParticles, Physics::m_p);

    double R11 = sqrt(gsl_stats_variance(&(dist.xDist_m[0]), 1, dist.xDist_m.size())) * 1e3;
    double R22 = sqrt(gsl_stats_variance(&(dist.pxDist_m[0]), 1, dist.pxDist_m.size()));
    double R33 = sqrt(gsl_stats_variance(&(dist.yDist_m[0]), 1, dist.yDist_m.size())) * 1e3;
    double R44 = sqrt(gsl_stats_variance(&(dist.pyDist_m[0]), 1, dist.pyDist_m.size()));
    double R55 = sqrt(gsl_stats_variance(&(dist.tOrZDist_m[0]), 1, dist.tOrZDist_m.size())) * 1e3;
    double R66 = sqrt(gsl_stats_variance(&(dist.pzDist_m[0]), 1, dist.pzDist_m.size()));

    double R21 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.pxDist_m[0]), 1, dist.xDist_m.size()) * 1e3 /
                  (expectedR11 * expectedR22));
    double R43 = (gsl_stats_covariance(&(dist.yDist_m[0]), 1, &(dist.pyDist_m[0]), 1, dist.yDist_m.size()) * 1e3 /
                  (expectedR33 * expectedR44));
    double R51 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.tOrZDist_m[0]), 1, dist.xDist_m.size()) * 1e6 /
                  (expectedR11 * expectedR55));
    double R52 = (gsl_stats_covariance(&(dist.pxDist_m[0]), 1, &(dist.tOrZDist_m[0]), 1, dist.pxDist_m.size()) * 1e3 /
                  (expectedR22 * expectedR55));
    double R61 = (gsl_stats_covariance(&(dist.xDist_m[0]), 1, &(dist.pzDist_m[0]), 1, dist.xDist_m.size()) * 1e3 /
                  (expectedR11 * expectedR66));
    double R62 = (gsl_stats_covariance(&(dist.pxDist_m[0]), 1, &(dist.pzDist_m[0]), 1, dist.pxDist_m.size()) /
                  (expectedR22 * expectedR66));

    EXPECT_NEAR(R11 / expectedR11, 1.0, 0.002);
    EXPECT_NEAR(R22 / expectedR22, 1.0, 0.002);
    EXPECT_NEAR(R33 / expectedR33, 1.0, 0.002);
    EXPECT_NEAR(R44 / expectedR44, 1.0, 0.002);
    EXPECT_NEAR(R55 / expectedR55, 1.0, 0.002);
    EXPECT_NEAR(R66 / expectedR66, 1.0, 0.002);

    EXPECT_NEAR(R21 / expectedR21, 1.0, 0.02);
    EXPECT_NEAR(R43 / expectedR43, 1.0, 0.02);
    EXPECT_NEAR(R51 / expectedR51, 1.0, 0.05);
    EXPECT_NEAR(R52 / expectedR52, 1.0, 0.02);
    EXPECT_NEAR(R61 / expectedR61, 1.0, 0.02);
    EXPECT_NEAR(R62 / expectedR62, 1.0, 0.02);
}