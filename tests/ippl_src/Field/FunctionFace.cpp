#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Field/BareField.h"
#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Index/Index.h"

#include <iostream>

constexpr double roundOffError = 1e-10;

double neg  (const double& x) { return -x; }
double twice(const double& x) { return 2.0*x; }
double zero (const double&)   { return 0.0; }
double one  (const double&)   { return 1.0; }

TEST(Field, FunctionFace)
{
    OpalTestUtilities::SilenceTest silencer;
    const unsigned Dim=2;
    Index I(5),J(5);
    FieldLayout<Dim> layout(I,J);
    Field<double,Dim> B(layout);

    // Set initial boundary conditions.
    BConds<double,Dim> bc;
    bc[0] = new FunctionFace<double,Dim>(neg,0);
    bc[1] = new FunctionFace<double,Dim>(twice,1);
    bc[2] = new FunctionFace<double,Dim>(zero,2);
    bc[3] = new FunctionFace<double,Dim>(one,3);

    // An array for testing.
    Field<double,Dim> A(layout,GuardCellSizes<Dim>(1),bc);

    A = 0.0;
    A[I][J] += I + J*10;
    std::cout << A << std::endl;
    EXPECT_NEAR(sum(A),550,roundOffError);
    B[I][J] = A[I-1][J-1];
    std::cout << B << std::endl;
    EXPECT_NEAR(sum(B),204,roundOffError);
    B[I][J] = A[I+1][J+1];
    std::cout << B << std::endl;
    EXPECT_NEAR(sum(B),677,roundOffError);
}