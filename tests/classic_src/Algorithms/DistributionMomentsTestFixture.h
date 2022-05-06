#include "Algorithms/DistributionMoments.h"

#include "gtest/gtest.h"

class OpalParticle;

class DistributionMomentsTest: public ::testing::Test {
protected:
    static void SetUpTestCase();

    static void addParticle(double x, double px, double y, double py, double z, double pz);

    static std::vector<OpalParticle> particles_m;
    static DistributionMoments distributionMoments_m;
};