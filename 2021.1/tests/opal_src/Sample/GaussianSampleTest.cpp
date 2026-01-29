#include "gtest/gtest.h"

#include "Sample/SampleGaussianSequence.h"
#include "opal_test_utilities/SilenceTest.h"

TEST(GaussSampleTest, ChainTest) {
    OpalTestUtilities::SilenceTest silencer;
    unsigned int nSample = 101;
    SampleGaussianSequence seq(-5, 5, 1, nSample);

    unsigned int j = 0;
    for (unsigned int i = 0; i * 2 < nSample - 1; ++ i) {
        seq.getNext(i);
        ++j;
    }
    double x = seq.getNext(j);

    EXPECT_NEAR(x, 0.0, 1e-8);
}