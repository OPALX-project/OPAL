#include <math.h>

#include "gtest/gtest.h"
#include "Classic/AbsBeamline/EndFieldModel/AsymmetricEnge.h"

double engeTest(double delta_s, double c, double lambda_end) {
    double value = 1/(1+exp(c*delta_s/lambda_end));
    return value;
}

double tryEnge(double s, double c1, double c2, double l1, double l2, double x1, double x2) {
    double out = -1.0;
    out += engeTest(-s-x1, c1, l1);
    out += engeTest(-x2+s, c2, l2);
    return out;
}

TEST(AsymmetricEngeTest, FunctionTest) {
    endfieldmodel::AsymmetricEnge enge({0.0, 1.0}, 10.0, 3.0, 
                                       {0.0, 4.0}, 11.0, 6.0);
    for (double s = -20.0; s < 20.0; s += 1.0) {
        double test = enge.function(s, 0);
        double ref = tryEnge(s, 1, 4, 3, 6, 10, 11);
        EXPECT_NEAR(test, ref, 1e-12) << s << " " << test << " " << ref;
    }

}

TEST(AsymmetricEngeTest, FunctionTest2) {
    endfieldmodel::AsymmetricEnge enge({0.0, 1.0}, 10.0, 3.0e-3, // c, x0, lambda
                                       {0.0, 4.0}, 11.0, 6.0e-3);
    double test1 = enge.function(11.0, 0);
    EXPECT_NEAR(test1, 0.5, 1e-3) << test1;
    double test2 = enge.function(-10.0, 0);
    EXPECT_NEAR(test2, 0.5, 1e-3) << test2;
}

TEST(AsymmetricEngeTest, DerivativeTest) {
    endfieldmodel::AsymmetricEnge enge({0.0, 1.0}, 10.0, 3.0, 
                                       {0.0, 4.0}, 11.0, 6.0);
    double ds = 0.01;
    for (int n = 1; n < 5; ++n) {
        for (double s = -20.0; s < 20.0; s += 1.0) {
            double fplus = enge.function(s+ds, n-1);
            double fminus = enge.function(s-ds, n-1);
            double test = enge.function(s, n);
            double ref = (fplus-fminus)/2/ds;
            double delta = abs(ref*1e-3);
            EXPECT_NEAR(test, ref, delta) << s << " " << test << " " << ref;
        }
    }
}
