#include <cmath>
#include "gtest/gtest.h"
#include "Classic/AbsBeamline/EndFieldModel/Enge.h"

TEST(EngeTest, DerivativeTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 10.0, 0.5);
    enge.setMaximumDerivative(0);
    EXPECT_NEAR(enge.function(10.0, 0), enge.function(-10.0, 0), 1e-12);
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        enge.setMaximumDerivative(i+1);
        double dfdxNumerical = (enge.function(10.0+dx, i)-
                                enge.function(10.0-dx, i))/2/dx;
        double dfdx = enge.function(10.0, i+1);
        EXPECT_NEAR(dfdx/dfdxNumerical, 1.0, 1e-8   ) << " for " << i << "^th derivative";
    }
}

double myEnge(double x, std::vector<double> a, double x0, double lambda) {
    double deltaX = (x-x0)/lambda;
    double xPow = 1.0;
    double p = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        p += a[i]*xPow;
        xPow *= deltaX;
    }
    double enge = 1/(1.0+std::exp(p));
    return enge;
}

TEST(EngeTest, FunctionTest) {
    std::vector<double> zVector = {0.0, 1.0, 2.0, 3.0};
    endfieldmodel::Enge enge({0.0, 1.0}, 0.0, 0.5);
    for (auto z: zVector) {
        EXPECT_NEAR(enge.getEnge(z, 0), myEnge(z, {0.0, 1.0}, 00.0, 0.5), 1e-12);
    }
}

TEST(EngeTest, HNTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 10.0, 0.5);
    enge.setMaximumDerivative(11);
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        double dhdxNumerical = (enge.hN(10.0+dx, i)-
                                enge.hN(10.0-dx, i))/2/dx;
        double dhdx = enge.hN(10.0, i+1);
        EXPECT_NEAR(dhdx, dhdxNumerical, 1e-5)
                << " for " << i << "^th derivative";
    }
}

TEST(EngeTest, GNTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 1.0, 0.5);
    enge.setMaximumDerivative(11);
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        double dgdxNumerical = (enge.gN(0.1+dx, i)-
                                enge.gN(0.1-dx, i))/2/dx;
        double dgdx = enge.gN(0.1, i+1);
        EXPECT_NEAR(dgdx/dgdxNumerical, 1.0, 1e-5)
                << " for " << i << "^th derivative";
    }
}

TEST(EngeTest, RescaleTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge1 = endfieldmodel::Enge(a, 10.0, 0.5);
    enge1.setMaximumDerivative(0);
    endfieldmodel::Enge* enge2 = enge1.clone();
    ASSERT_NE(enge2, nullptr);
    EXPECT_EQ(enge1.function(10.0, 0), enge2->function(10.0, 0));
    enge2->rescale(0.1);
    EXPECT_EQ(enge1.function(10.0, 0), enge2->function(1.0, 0));
    double diffScale = 1.0;
    for(size_t i = 0; i < 10; ++i) {
        double dfdx1 = enge1.function(10.0, i);
        double dfdx2 = enge2->function(1.0, i);
        EXPECT_NEAR(dfdx1*diffScale/dfdx2, 1.0, 1e-12)
                    << " for " << i << "^th derivative";
        diffScale *= 10.0;
    }
    delete enge2;
}

