#include "gtest/gtest.h"

#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/CartesianCentering.h"
#include "AppTypes/Vektor.h"

#include "opal_test_utilities/SilenceTest.h"

#include <cstdio>
#include <iostream>

constexpr double roundOffError = 1e-10;

// template definition in BCond.cpp
//CenteringEnum CCCEnums<2U,2U,0U>::vectorFace[2U*2U];

TEST(Field, Eureka)
{
    OpalTestUtilities::SilenceTest silencer;

    const int N = 6;
    Index I(1,N), J(1,N);
    Index I0(-1,N+2), J0(-1,N+2);
    Index I1(2,N-1),J1(2,N-1);

    FieldLayout<2> layout(I,J);
    FieldLayout<2> layout0(I0,J0);
    GuardCellSizes<2> gc(2);

    // Test all cell centering.
    BConds<double,2> bc1;
    bc1[1] = new EurekaFace<double,2>(0);
    bc1[2] = new EurekaFace<double,2>(1);
    bc1[3] = new EurekaFace<double,2>(2);
    bc1[4] = new EurekaFace<double,2>(3);

    Field<double,2> A1(layout,gc,bc1), A0(layout0);

    A1[I][J] = 10.0*I + 100.0*J;

    // Assign the full domain to A0.
    A0[I0][J0] = A1[I0][J0];

    // See if we got the right answer.
    // s1 should be zero.
    A0[I1][J1] -= 10.0*I1 + 100.0*J1;

    double s1 = sum(A0*A0);

    EXPECT_NEAR(s1, 0.0, roundOffError);

    // Test mixed centering.
    typedef Vektor<double,2> T;
    typedef UniformCartesian<2> M;
    typedef CartesianCentering<CCCEnums<2,2,0>::vectorFace,2,2> C;
    BConds<T,2,M,C> bc2;
    bc2[0] = new EurekaFace<T,2,M,C>(0,0);
    bc2[1] = new EurekaFace<T,2,M,C>(1,0);
    bc2[2] = new EurekaFace<T,2,M,C>(2,0);
    bc2[3] = new EurekaFace<T,2,M,C>(3,0);
    bc2[4] = new EurekaFace<T,2,M,C>(0,1);
    bc2[5] = new EurekaFace<T,2,M,C>(1,1);
    bc2[6] = new EurekaFace<T,2,M,C>(2,1);
    bc2[7] = new EurekaFace<T,2,M,C>(3,1);

    Field<T,2,M,C> B1(layout,gc,bc2), B0(layout0);

    // Fill with nontrivial data.
    B1[I][J] = T(1,1)*(I + 10.0*J);

    std::cout << B1 << std::endl << std::endl;

    // Pull it out into a field that shows the guard layers.
    B0[I0][J0] = B1[I0][J0];

    std::cout << B0 << std::endl << std::endl;

    // See if we got the right answer.
    B0[I1][J1] -= T(1,1)*(I1+10.0*J1);
    B0[1][J1]  -= T(1,1)*(1+10.0*J1);
    B0[N][J1]  -= T(1,1)*(N+10.0*J1);
    B0[I1][1]  -= T(1,1)*(I1+10.0);
    B0[I1][N]  -= T(1,1)*(I1+10.0*N);

    B0[N][N]  -= T(1,1)*(N+10.0*N);
    B0[1][N]  -= T(1,1)*(1+10.0*N);
    B0[1][1]  -= T(1,1)*(1+10.0);
    B0[N][1]  -= T(1,1)*(N+10.0);

    Vektor<T,2> s2 = sum(B0*B0);
    EXPECT_TRUE((s2 == Vektor<T,2>(0,0)) );

    std::cout << s2 << std::endl;
    std::cout << B0 << std::endl << std::endl;
}