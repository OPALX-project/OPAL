#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Index/SOffset.h"

#include <iostream>

constexpr unsigned Dim = 2;
constexpr double   roundOffError = 1e-10;

TEST(Index, SOffset)
{
    // testing SOffset creation
    SOffset<Dim> A(1, 2);
    SOffset<Dim> B(1, 1);
    EXPECT_NEAR(A[0],1, roundOffError);
    EXPECT_NEAR(A[1],2, roundOffError);
    EXPECT_NEAR(B[0],1, roundOffError);
    EXPECT_NEAR(B[1],1, roundOffError);

    // testing SOffset [] operator
    SOffset<Dim> soLeft;
    SOffset<Dim> soRight;
    for (unsigned int d=0; d < Dim; d++) {
        soLeft[d] = 0;
        soRight[d] = (d == 1 ? 1 : 0);
    }
    EXPECT_NEAR(soLeft[0], 0, roundOffError);
    EXPECT_NEAR(soLeft[1], 0, roundOffError);
    EXPECT_NEAR(soRight[0], 0, roundOffError);
    EXPECT_NEAR(soRight[1], 1, roundOffError);

    // testing SOffset addition, subtraction, and copy constructor
    SOffset<Dim> C(A + B);
    SOffset<Dim> D(A - B);
    EXPECT_NEAR(C[0], 2, roundOffError);
    EXPECT_NEAR(C[1], 3, roundOffError);
    EXPECT_NEAR(D[0], 0, roundOffError);
    EXPECT_NEAR(D[1], 1, roundOffError);

    // testing SOffset +=, -= operators
    D += C;
    EXPECT_NEAR(D[0], 2, roundOffError);
    EXPECT_NEAR(D[1], 4, roundOffError);
    C -= D;
    EXPECT_NEAR(C[0],  0, roundOffError);
    EXPECT_NEAR(C[1], -1, roundOffError);

    // testing SOffset comparisons
    EXPECT_TRUE(B < A);
    EXPECT_TRUE(B <= A);
    EXPECT_TRUE(D > C);
    EXPECT_TRUE(D >= C);
    EXPECT_TRUE(A != B);
    EXPECT_TRUE(B == B);

    // testing containment check
    NDIndex<Dim> N1(Index(0,1), Index(0,1));
    EXPECT_TRUE(! A.inside(N1));
    EXPECT_TRUE(  B.inside(N1));
}

#include "Index/SIndex.h"

TEST(Index, SIndex)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(25);
    Index J(50);
    Index I2(2);
    Index J2(5);
    NDIndex<Dim> NDX(I,J);
    int IP[Dim];
    IP[0] = -1;
    IP[1] = -1;

    SOffset<Dim> so1(1,1);
    SOffset<Dim> so2(IP);
    std::cout << "Created SOffset so1 = " << so1 << std::endl;
    std::cout << "Created SOffset so2 = " << so2 << std::endl;
    std::cout << "Adding IP to so1  = " << so1 + IP << std::endl;
    EXPECT_EQ((so1+IP)[0],0);
    std::cout << "Adding so1 to IP = " << IP + so1 << std::endl;
    EXPECT_EQ((IP+so1)[1],0);
    std::cout << "Adding so2 to so1 = " << so1 + so2 << std::endl;
    EXPECT_EQ((so1+so2)[0],0);
    std::cout << "Adding so1 to NDX " << NDX << " = " << NDX + so1 << std::endl;
    std::cout << "Adding NDX " << NDX << " to so1 = " << so1 + NDX << std::endl;
    std::cout << "Multiplying NDX " << NDX << " by so2 = " << NDX * so2 << std::endl;
    std::cout << "Multiplying so2 by NDX " << NDX << " = " << so2 * NDX  << std::endl;
    so1 -= IP;
    so2 += IP;
    std::cout << "Accumulated IP from so1 = " << so1 << std::endl;
    EXPECT_EQ((so1)[0],2);
    std::cout << "Accumulated IP into so2 = " << so2 << std::endl;
    EXPECT_EQ((so2)[0],-2);

    FieldLayout<Dim> layout(I, J, PARALLEL, PARALLEL, 2*Ippl::getNodes());
    Field<double,Dim> A(layout);
    Field<bool,Dim> B(layout);

    SIndex<Dim> s1(layout);
    SIndex<Dim> s2 = s1(1,-1);
    SIndex<Dim> x3 = s1(IP);

    std::cout << "Created s1 = " << s1 << std::endl;
    std::cout << "Created s2 = " << s2 << std::endl;
    std::cout << "Created x3 = " << x3 << std::endl;

    s1.addIndex(NDIndex<Dim>(Index(2), Index(3)));
    s2.addIndex(SOffset<Dim>(0,0));
    s2.addIndex(NDIndex<Dim>(Index(20,23), Index(45,46)));

    // testmsg << "Added new points, s1 = "; s1.printDebug(testmsg);
    // testmsg << "Added new points, s2 = "; s2.printDebug(testmsg);

    SIndex<Dim> s3(s1);
    s3 = NDIndex<Dim>(Index(1,5), Index(2,4));

    std::cout << "Created s3 = " << s3 << std::endl;

    s3 &= s1;

    std::cout << "Intersection of s3 and s1 = " << s3 << std::endl;

    // now, test assigment of a Field expression to an SIndex
    A[I][J] = I + J;
    B = lt(A,10);
    s3 = lt(A,10);

    // do a union with a slightly different condition
    s3 |= (gt(A,15) && lt(A,20));
    std::cout << "union of s3 and where 15 < A < 20 = " << s3 << std::endl;

    // do an intersection with an overlapping condition
    s3 &= (gt(A,8) && lt(A,12));
    std::cout << "intersection of s3 and where 8 < A < 12 = " << s3 << std::endl;

    // do an indexed assignment to s3
    s3[I2][J2] = (lt(A[I2 + 1][J2+1], 5) && gt(A[I2 + 2][J2 + 2], 0));
    std::cout << "s3[I2][J2] = expr ==> s3 = " << s3 << std::endl;

    // now use s3 in a Field expression
    std::cout << "Originally, A = I + J ... after A[s3] = A[s3(1,1)]: A = ";
    A[s3] = A[s3(1,1)] + 10;
    std::cout << A << std::endl;

    EXPECT_EQ(sum(A), 45685);
}

TEST(Index, SubField)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(4);
    Index I2(2);
    Index J(4);
    Index K(2);
    FieldLayout<Dim> layout(I, J, PARALLEL, PARALLEL, 4*Ippl::getNodes());
    FieldLayout<3> layout3(I, J, K, PARALLEL, PARALLEL, PARALLEL,
                           2*Ippl::getNodes());
    BConds<double,Dim> dbc;
    for (unsigned int f=0; f < 2*Dim; f++) dbc[f] = new ZeroFace<double,Dim>(f);
    Field<double,Dim> A(layout, dbc, GuardCellSizes<Dim>(1));
    A = 0.0;
    BConds<bool,Dim> bbc;
    for (unsigned int f=0; f < 2*Dim; f++) bbc[f] = new ZeroFace<bool,Dim>(f);
    Field<bool,Dim>   B(layout, bbc, GuardCellSizes<Dim>(1));
    B = true;
    Field<double,Dim> C(layout, dbc, GuardCellSizes<Dim>(1));
    C = 0.0;

    BConds<float,3> fbc;
    for (unsigned int f=0; f < 2*3; f++) fbc[f] = new ZeroFace<float,3>(f);
    Field<float,3>    A3(layout3, fbc, GuardCellSizes<3>(2));
    A3 = 0.0;
    Field<float,3>    B3(layout3, fbc, GuardCellSizes<3>(2));
    B3 = 0.0;
    SIndex<Dim> s1(layout);
    s1.addIndex(NDIndex<Dim>(Index(1,1), Index(1,2)));

    SIndex<Dim> s2 = s1(1,-1);
    s2.addIndex(SOffset<Dim>(0,0));

    std::cout << "\n************ testing 2D SubParticleAttrib<> ************" << std::endl;
    A[I][J] = (I+1) + (J*10);
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << " A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),280);
    // not working anymore
    // ParticleAttrib<double> PA;
    // PA[s1] = A[s1];
    // std::cout << "Result of PA[s1] = A[s1] : PA = " << std::endl;
    // std::cout << PA[s1] << std::endl;

    // std::cout << "s2 = " << s2 << std::endl;
    // PA[s2] = A[s2] + A[s2];
    // std::cout << "Result of PA[s2] = A[s2] + A[s2] : PA = " << std::endl;
    // std::cout << PA[s2] << std::endl;

    // PA[s2] *= (A[s2] + PA[s2]);
    // std::cout << "Result of PA[s2] *= A[s2] + PA[s2] : PA = " << std::endl;
    // std::cout << PA[s2] << std::endl;

    // A = 1.0;
    SIndex<Dim> sX = s1;
    sX.addIndex(SOffset<Dim>(0,0));
    std::cout << "Created sX = " << sX << std::endl;

    // A[sX] += PA[sX];
    // std::cout << "Result of A[sX] = 1 + PA[sX] : A = " << A << std::endl;

    A = 0.0;
    C = 1.0;
    A[sX] = C[sX];
    std::cout << "Result of A[sX] = C[sX], after C = 1.0: A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),3);

    // ParticleAttrib<Vektor<double,2> > PA2;
    // PA2(0)[sX] = A[sX] + PA[sX];
    // PA2(1)[sX] = -(A[sX] + PA[sX]);
    // std::cout << "Result of PA2[sX] = +- (A[sX] + PA[sX]):" << std::endl;
    // std::cout << "  A = " << A << std::endl;
    // std::cout << " PA = " << PA[sX] << std::endl;
    // std::cout << "PA2 = " << PA2[sX] << std::endl;

    std::cout << "\n************ testing 2D SubBareField<NDIndex> ************" << std::endl;

    Index IX(0,2);
    Index JX(2,2);
    NDIndex<Dim> sub1(IX, JX);
    SubBareField<double,Dim,NDIndex<Dim> > SA1(A, sub1);
    std::cout << "Created SubBareField SA1 from A and NDIndex " << sub1 << ":" <<std::endl;
    std::cout << " SA1 = " << SA1 << std::endl;

    SA1 = 3.0;
    std::cout << "SA1 set to 3:" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),11);

    NDIndex<Dim> sub2(IX + 1, JX);
    SubBareField<double,Dim,NDIndex<Dim> > SA2(A, sub2);
    std::cout << "Created SubBareField SA2 from A and NDIndex " << sub2 << ":" <<std::endl;
    std::cout << " SA2 = " << SA2 << std::endl;

    SA1 = SA2 + 4.0;
    std::cout << "SA1 on " << sub1 << " set to (SA2 + 4):" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),20);

    std::cout << "\n************ testing 2D SubBareField<SIndex> ************" << std::endl;

    SubBareField<double,Dim,SIndex<Dim> > SB1(A, s1);
    std::cout << "Created SubBareField SB1 from A and SIndex " << s1 << ":" << std::endl;
    std::cout << " SB1 = " << SB1 << std::endl;

    SB1 = -2.0;
    std::cout << "SB1 set to -2:" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),5);

    SubBareField<double,Dim,SIndex<Dim> > SB2(A, s1(1,-1));
    std::cout << "Created SubBareField SB2 from A and SIndex " << SB2.getDomain();
    std::cout << ":" << std::endl;
    std::cout << " SB2 = " << SB2 << std::endl;

    SB1 = SB2 + 4.0;
    std::cout << "SB1 set to (SB2 + 4):" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),23);

    std::cout << "\n************ testing 2D SubBareField<SOffset> ************" << std::endl;

    SubBareField<double,Dim,SOffset<Dim> > SC1(A, SOffset<Dim>(2,1));
    std::cout << "Created SubBareField SC1 from A and SOffset "  << SC1.getDomain();
    std::cout << ":" << std::endl;
    std::cout << " SC1 = " << SC1 << std::endl;

    SC1 = 10.0;
    std::cout << "SC1 set to 10:" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),33);

    SubBareField<double,Dim,SOffset<Dim> > SC2(A, SOffset<Dim>(2,2));
    SubBareField<double,Dim,SOffset<Dim> > SC3(A, SOffset<Dim>(3,2));
    SC1 = SC1 + (SC2 + SC3) * 10;
    std::cout << "SC1 set to 10 + [(2,2) + (3,2)]*10:" << std::endl;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),73);

    std::cout << "\n************ testing 2D Field[SIndex] ************" << std::endl;

    std::cout << "Current SIndex s1 = " << s1 << std::endl;
    std::cout << "Results of A[s1] = -1000:" << std::endl;
    A[s1] = -1000.0;
    std::cout << "   A = " << A << " " << sum(A) << std::endl;
    EXPECT_EQ(sum(A),-2939);
    B[s1] = ne(A[s1], A[s1]);
    std::cout << "Results of B[s1] = (A[s1] != A[s1]):" << std::endl;
    std::cout << "   B = " << B << " " << sum(B) << std::endl;
    EXPECT_EQ(sum(B),1);
    C[s1] = 1.0;
    std::cout << "Results of C[s1] = 1:" << std::endl;
    std::cout << "   C = " << C << " " << sum(C) << std::endl;
    EXPECT_EQ(sum(C),16);
    C[s1] = A[s1(-1,-1)] - 2000;
    std::cout << "Results of C[s1] = (A[s1(-1,-1)] - 2000) :" << std::endl;
    std::cout << "   C = " << C << " " << sum(C) << std::endl;
    EXPECT_EQ(sum(C),-6987);

    std::cout << "\n************ testing 2D Field[SIndex,compressed] ****" << std::endl;

    Field<double,Dim> A2(layout);
    A2 = 1.0;
    std::cout << "   A2 = " << A2 << " " << sum(A2) << std::endl;
    EXPECT_EQ(sum(A2),16);
    std::cout << "   A2.compressedFraction = " << A2.CompressedFraction() << std::endl;
    A2[I2][I2] = -10.0;
    std::cout << "Initial settings:" << std::endl;
    std::cout << "   A2 = " << A2 << " " << sum(A2) << std::endl;
    EXPECT_EQ(sum(A2),-28);
    std::cout << "   A2.compressedFraction = " << A2.CompressedFraction() << std::endl;
    s1 = ne(A2, 0.0);
    std::cout << "   s1 = " << s1 << std::endl;
    std::cout << "   A2 = " << A2 << " " << sum(A2) << std::endl;
    EXPECT_EQ(sum(A2),-28);
    std::cout << "   A2.compressedFraction = " << A2.CompressedFraction() << std::endl;
    A2[s1] = A2[s1] + 100.0;
    std::cout << "Results of A2[s1] = A2[s1] + 100.0 :" << std::endl;
    std::cout << "   A2 = " << A2 << " " << sum(A2) << std::endl;
    EXPECT_EQ(sum(A2),1572);
    std::cout << "   A2.compressedFraction = " << A2.CompressedFraction() << std::endl;
    std::cout << "   s1 = " << s1 << std::endl;

    std::cout << "\n************ testing 3D Field[SIndex] ************" << std::endl;

    A3[I][J][K] = I + J + K;
    B3 = 0.0;

    SIndex<3> sindx3(layout3);
    sindx3 = eq(A3, 2.0);

    SOffset<3> offset;
    offset[0] = 1;
    offset[1] = 0;
    offset[2] = (-1);

    std::cout << "Current 3D SIndex sindx3 = " << sindx3 << std::endl;
    std::cout << "Current 3D A3 = " << A3 << " " << sum(A3) << std::endl;
    EXPECT_EQ(sum(A3),112);
    std::cout << "A3 after A3[sindx3] = -1:" << std::endl;
    A3[sindx3] = -1.0;
    std::cout << "   A3 = " << A3 << " " << sum(A3) << std::endl;
    EXPECT_EQ(sum(A3),97);
    B3[sindx3] = 1.0;
    std::cout << "Results of B3[sindx3] = 1:" << std::endl;
    std::cout << "   B3 = " << B3 << " " << sum(B3) << std::endl;
    EXPECT_EQ(sum(B3),5);
    B3[sindx3] = A3[sindx3(offset)] - 2000.0;
    std::cout << "Results of B3[sindx3] = A3[sindx3(1,0,-1)] - 2000.0:" << std::endl;
    std::cout << "   B3 = " << B3 << " " << sum(B3) << std::endl;
    EXPECT_EQ(sum(B3),-10002);
}

#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"

/***************************************************************************
  A simple program to test SubField assignments using two fields, one on
  a cell-centered layout, the other on a vertex-centered (but otherwise
  similar) layout.
***************************************************************************/

TEST(Index, VertCell)
{
    OpalTestUtilities::SilenceTest silencer;


    Index I(4), J(4), U(5), V(5);
    FieldLayout<Dim> LC(I, J, PARALLEL, PARALLEL, 4);
    FieldLayout<Dim> LV(U, V, PARALLEL, PARALLEL, 4);

    Field<int,Dim,UniformCartesian<Dim>,Cell> A(LC, GuardCellSizes<Dim>(1));
    Field<int,Dim,UniformCartesian<Dim>,Vert> B(LV, GuardCellSizes<Dim>(1));

    A[I][J] = J + 1;
    B[U][V] = -V - 1;
    std::cout << "A at start: " << A << std::endl;
    std::cout << "B at start: " << B << std::endl;

    SIndex<Dim> SI(LC);
    SI = eq(A,2) || eq(A,4);
    std::cout << "SI at start: " << SI << std::endl;

    A[SI] = A[SI(0,-1)] + B[SI(1,1)];

    std::cout << "A at end: " << A << std::endl;

    int s = sum(A);

    EXPECT_EQ(s, 0);
}