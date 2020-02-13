#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Field/Field.h"
#include "Field/BCond.h"
#include "Utility/FieldDebug.h"
#include "Index/Index.h"
#include "Particle/ParticleSpatialLayout.h"

#include <iostream>

constexpr double roundOffError = 1e-10;

typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;

TEST(Field, PeriodicBC)
{
    OpalTestUtilities::SilenceTest silencer;

    constexpr unsigned Dim = 3;

    Index I(3);
    Index J(3);
    Index K(3);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    domain[2] = K;
    FieldLayout<Dim> layout(domain);
    typedef UniformCartesian<Dim> M;

    // Set Cell-centered boundary conditions.
    BConds<double,Dim,M,Cell> cbc;
    cbc[0] = new ZeroFace<double,Dim,M,Cell>(0);
    cbc[1] = new ZeroFace<double,Dim,M,Cell>(1);
    cbc[2] = new ZeroFace<double,Dim,M,Cell>(2);
    cbc[3] = new ZeroFace<double,Dim,M,Cell>(3);
    cbc[4] = new ParallelPeriodicFace<double,Dim,M,Cell>(4);
    cbc[5] = new ParallelPeriodicFace<double,Dim,M,Cell>(5);

    BConds<Vector_t,Dim,M,Cell> vcbc;
    vcbc[0] = new ZeroFace<Vector_t,Dim,M,Cell>(0);
    vcbc[1] = new ZeroFace<Vector_t,Dim,M,Cell>(1);
    vcbc[2] = new ZeroFace<Vector_t,Dim,M,Cell>(2);
    vcbc[3] = new ZeroFace<Vector_t,Dim,M,Cell>(3);
    vcbc[4] = new ParallelPeriodicFace<Vector_t,Dim,M,Cell>(4);
    vcbc[5] = new ParallelPeriodicFace<Vector_t,Dim,M,Cell>(5);

    std::cout << "++++++++BConds object cbc begin++++++++" << std::endl;
    std::cout << cbc;
    std::cout << "++++++++BConds object cbc end++++++++++" << std::endl;

    // Cell-centered Field's:
    std::cout << "layout: " << layout << std::endl;
    Field<double,Dim,M,Cell> cA(layout,GuardCellSizes<Dim>(1),cbc);
    Field<double,Dim,M,Cell> cB(layout,GuardCellSizes<Dim>(1),cbc);

    Field<Vector_t,Dim,M,Cell> dVf(layout,vcbc,GuardCellSizes<Dim>(1));

    // Assign reference values:
    int i,j,k;
    for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
            for(k=0; k<3; k++) {
                if (i==1 && j==1 && k==1)
                    assign(cA[i][j][k], 1.0);
                else
                    assign(cA[i][j][k], 0.0);
            }
        }
    }
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, Cell-centered case:

    std::cout << "++++++++++cA+++++++++++" << std::endl ;
    fp3(cA);
    EXPECT_NEAR(sum(cA),1.0,roundOffError);

    cB[I][J][K] = cA[I][J][K-2];
    std::cout << "++++++++++cB+++++++++++" << std::endl ;
    fp3(cB);
    EXPECT_NEAR(sum(cB),0.0,roundOffError);

    dVf = Grad(cA,dVf);
    std::cout << "++++++++++Grad(cA)+++++" << std::endl ;
    fp3(dVf);
    std::cout << sum(dVf) << std::endl;
    EXPECT_NEAR(sum(dVf)[0],0,roundOffError);
    EXPECT_NEAR(sum(dVf)[1],0,roundOffError);
    EXPECT_NEAR(sum(dVf)[2],0,roundOffError);

    dVf = Vector_t(0.0);

    dVf = Grad(cB,dVf);
    std::cout << "++++++++++Grad(cB)+++++" << std::endl ;
    fp3(dVf);
    std::cout << sum(dVf) << std::endl;
    EXPECT_NEAR(sum(dVf)[0],0,roundOffError);
    EXPECT_NEAR(sum(dVf)[1],0,roundOffError);
    EXPECT_NEAR(sum(dVf)[2],0,roundOffError);
}
