#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Field/BareField.h"
#include "Field/Field.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "FieldLayout/FieldLayout.h"
#include "Index/Index.h"
#include "Meshes/UniformCartesian.h"

#include <map>
#include <iostream>
#include <vector>

constexpr unsigned Dim = 2;
constexpr unsigned D1 = 1;
constexpr unsigned D2 = 2;
constexpr unsigned D3 = 3;
constexpr unsigned D4 = 4;
constexpr double   roundOffError = 1e-10;

TEST(Field, Balance)
{
    OpalTestUtilities::SilenceTest silencer;

    const int N=10;
    Index I(N);
    Index J(N);
    FieldLayout<Dim> layout1(I,J,PARALLEL,SERIAL,4);
    FieldLayout<Dim> layout2(I,J,SERIAL,PARALLEL,8);
    Field<int,Dim> A1(layout1),A2(layout2);

    A1 = 0;
    A1[I][J] += I + 10*J;
    A2 = A1;

    A2[I][J] -= I + 10*J;
    A2 *= A2;
    int s = sum(A2);

    EXPECT_NEAR(s, 0, roundOffError);
}

TEST(Field, Bool)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(10);

    double amin;
    FieldLayout<D1>  layout(I);
    Field<double,D1,UniformCartesian<D1>,
        UniformCartesian<D1>::DefaultCentering> A(layout);
    Field<bool,D1,UniformCartesian<D1>,
      UniformCartesian<D1>::DefaultCentering>   B(layout);
    Field<double,D1,UniformCartesian<D1>,
        UniformCartesian<D1>::DefaultCentering> C(layout);

    A[I] = I;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),45,roundOffError);

    amin = 5.0;
    C = amin;
    std::cout << "sum C = " << sum(C) << std::endl;
    EXPECT_NEAR(sum(C),50,roundOffError);

    B = lt( A, amin ) ;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),1,roundOffError);

    A = where( B, C, A ) ;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),60,roundOffError);

    A = where( lt(A,amin), A, C ) ;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),50,roundOffError);
}

TEST(Field, Component)
{
    Index I(5),J(5);
    FieldLayout<Dim> layout(I,J);
    typedef UniformCartesian<Dim,double> Mesh;
    Field<double,Dim,Mesh> S1(layout),S2(layout),S3(layout);
    Field<Vektor<double,Dim>,Dim,Mesh> V1(layout);

    S1[I][J] = I+10*J ;
    S2[I][J] = -I-10*J ;
    V1[I][J](0) << S1[I][J] ;
    V1[I][J](1) << S2[I][J]*10.0 ;
    S1[I][J] = V1[I][J](0) - (I+10*J);
    S2[I][J] = V1[I][J](1)/10.0 + (I+10*J);
    S1 *= S1;
    S2 *= S2;
    double s1 = sum(S1);
    double s2 = sum(S2);
    EXPECT_NEAR(s1,0,roundOffError);
    EXPECT_NEAR(s2,0,roundOffError);

    V1[I][J](0)  << I ;
    V1[I][J](1)  << 27.5 ;
    S1[I][J] = V1[I][J](0) - I;
    S2[I][J] = V1[I][J](1) - 27.5;
    S1 *= S1;
    S2 *= S2;
    s1 = sum(S1);
    s2 = sum(S2);
    EXPECT_NEAR(s1,0,roundOffError);
    EXPECT_NEAR(s2,0,roundOffError);

    S1[I][J] = I+10*J ;
    S2[I][J] = -I-10*J ;
    V1 = Vektor<double,2>(0,0);
    V1[I][J](0) += S1[I][J] ;
    V1[I][J](1) += S2[I][J]*10.0 ;
    S1[I][J] = V1[I][J](0) - (I+10*J);
    S2[I][J] = V1[I][J](1)/10.0 + (I+10*J);
    S1 *= S1;
    S2 *= S2;
    s1 = sum(S1);
    s2 = sum(S2);
    EXPECT_NEAR(s1,0,roundOffError);
    EXPECT_NEAR(s2,0,roundOffError);

    V1 = Vektor<double,2>(0,0);
    V1[I][J](0) += I ;
    V1[I][J](1) += 27.5 ;
    S1[I][J] = V1[I][J](0) - I;
    S2[I][J] = V1[I][J](1) + (-27.5);
    S1 *= S1;
    S2 *= S2;
    s1 = sum(S1);
    s2 = sum(S2);
    EXPECT_NEAR(s1,0,roundOffError);
    EXPECT_NEAR(s2,0,roundOffError);
}

TEST(Field, Compressed)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(10),J(10);
    FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL,4);
    GuardCellSizes<Dim> gc(1);
    typedef BareField<int,Dim> F;
    F A(layout,gc);

    // A should be constructed compressed.
    F::iterator_if lf;
    int count;

    //////////////////////////////////////////////////////////////////////
    // Test if it is constructed compressed
    // (or uncompressed, if --nofieldcompression)
    //////////////////////////////////////////////////////////////////////

    //    if (IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( ( (*lf).second->IsCompressed() &&  IpplInfo::noFieldCompression) ||
             (!(*lf).second->IsCompressed() && !IpplInfo::noFieldCompression)) {
            std::cout << "FAILED: An LField is (un)compressed," << count << std::endl;
            EXPECT_TRUE(false);
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Test whether fillGuardCells destroys the compression.
    //////////////////////////////////////////////////////////////////////
    A.fillGuardCells();
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( ( (*lf).second->IsCompressed() &&  IpplInfo::noFieldCompression) ||
             (!(*lf).second->IsCompressed() && !IpplInfo::noFieldCompression)) {
            std::cout << "FAILED: (un)compressed after fillGuardCells, " << count << std::endl;
            EXPECT_TRUE(false);
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning an index uncompresses.
    //////////////////////////////////////////////////////////////////////
    assign(A[I][J] , I + 10*J);
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( (*lf).second->IsCompressed() ) {
            std::cout << "FAILED: Compressed after assigning Index, " << count << std::endl;
            EXPECT_TRUE(false);
        }
    }

    int s = sum(A);
    int il = I.length();
    int jl = J.length();
    int ss = jl*il*(il-1)/2 + il*jl*(jl-1)*5;
    EXPECT_NEAR(s, ss, roundOffError);

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning a constant compresses.
    //////////////////////////////////////////////////////////////////////
    A = 1;
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( ( (*lf).second->IsCompressed() &&  IpplInfo::noFieldCompression) ||
             (!(*lf).second->IsCompressed() && !IpplInfo::noFieldCompression)) {
            std::cout << "FAILED: (Un)compressed after assigning constant, " << count << std::endl;
            EXPECT_TRUE(false);
        }
    }
    s = sum(A);
    EXPECT_EQ(s, (int)(I.length()*J.length()));

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning a constant with indexes compresses.
    //////////////////////////////////////////////////////////////////////
    // First uncompress it.
    assign(A[I][J] , I+10*J);
    // Then compress it.
    assign(A[I][J] , 1 );
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( ( (*lf).second->IsCompressed() &&  IpplInfo::noFieldCompression) ||
             (!(*lf).second->IsCompressed() && !IpplInfo::noFieldCompression)) {
            std::cout << "FAILED: (Un)compressed after I-assigning constant, " << count << std::endl;
            EXPECT_TRUE(false);
        }
    }
    s = sum(A);
    EXPECT_EQ(s, (int)(I.length()*J.length()));

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning a subrange uncompresses.
    //////////////////////////////////////////////////////////////////////
    Index I1( I.min()+1 , I.max()-1 );
    Index J1( J.min()+1 , J.max()-1 );
    A=0;
    assign(A[I1][J1] , 1 );
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
        if ( (*lf).second->IsCompressed() ) {
            std::cout << "FAILED: Compressed after assigning subrange, " << count << std::endl;
            EXPECT_TRUE(false);
        }
    }
    s = sum(A);
    EXPECT_EQ(s, (int)(I1.length()*J1.length()));

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning a single element uncompresses one.
    //////////////////////////////////////////////////////////////////////
    A = 0;
    if (!IpplInfo::noFieldCompression) {
        for (lf=A.begin_if(); lf!=A.end_if(); ++lf) {
            if( !(*lf).second->IsCompressed() ) {
                EXPECT_TRUE(false);
            }
        }
    }
    assign(A[3][3] , 1 );
    count = 0;
    if (IpplInfo::noFieldCompression) {
        if (A.CompressedFraction() != 0) {
            std::cout << "FAILED: Compressed somewhere after assigning single" << std::endl;
            EXPECT_TRUE(false);
        }
    } else {
        for (lf=A.begin_if(); lf!=A.end_if(); ++lf) {
            if ( (*lf).second->IsCompressed() )
                ++count;
        }
        int reduced_count = 0;
        reduce(&count,&count+1,&reduced_count,OpAddAssign());
        EXPECT_EQ(reduced_count, 3);
    }
    s = sum(A);
    EXPECT_NEAR(s, 1, roundOffError);

    //////////////////////////////////////////////////////////////////////
    // Test whether assigning a single element uncompresses one.
    //////////////////////////////////////////////////////////////////////
    A = 0;
    if (!IpplInfo::noFieldCompression) {
        for (lf=A.begin_if(); lf!=A.end_if(); ++lf) {
            if (! (*lf).second->IsCompressed() ) {
                EXPECT_TRUE(false);
            }
        }
    }
    assign(A[4][3] , 1 );
    count = 0;
    if (IpplInfo::noFieldCompression) {
        if (A.CompressedFraction() != 0) {
            std::cout << "FAILED: Compressed somewhere after assigning single" << std::endl;
            EXPECT_TRUE(false);
        }
    } else {
        for (lf=A.begin_if(); lf!=A.end_if(); ++lf) {
            if ( (*lf).second->IsCompressed() )
                ++count;
        }
        int reduced_count = 0;
        reduce(&count,&count+1,&reduced_count,OpAddAssign());
        if (Ippl::deferGuardCellFills) {
            EXPECT_EQ(reduced_count,3);
            std::cout << "PASSED: Uncompressed one after assigning in guard cell" << std::endl;
        }
        else {
            EXPECT_EQ(reduced_count,2);
            std::cout << "PASSED: Uncompressed two after assigning in guard cell" << std::endl;
        }
    }
    s = sum(A);
    EXPECT_NEAR(s, 1, roundOffError);

    //////////////////////////////////////////////////////////////////////
    // Test whether an operation on that array leaves it correct.
    //////////////////////////////////////////////////////////////////////
    A *= A;
    count = 0;
    if (IpplInfo::noFieldCompression) {
        EXPECT_EQ(A.CompressedFraction(),0);
        std::cout << "FAILED: Compressed somewhere after squaring" << std::endl;
    }
    else {
        for (lf=A.begin_if(); lf!=A.end_if(); ++lf) {
            if ( (*lf).second->IsCompressed() )
                ++count;
        }
        int reduced_count = 0;
        reduce(&count,&count+1,&reduced_count,OpAddAssign());
        if (Ippl::deferGuardCellFills) {
            EXPECT_EQ(reduced_count,3);
            std::cout << "PASSED: Uncompressed one after squaring" << std::endl;
        }
        else {
            EXPECT_EQ(reduced_count,2);
            std::cout << "PASSED: Uncompressed two after squaring" << std::endl;
        }
    }
    s = sum(A);
    EXPECT_NEAR(s, 1, roundOffError);

    //////////////////////////////////////////////////////////////////////
    // Make sure we can construct a Field of Maps.
    //////////////////////////////////////////////////////////////////////
    BareField< std::map<int,double> , Dim > B(layout);
    BareField< std::map<int,double> , Dim >::iterator_if lb;
    std::map<int,double> m;
    m[1] = cos(1.0);
    m[2] = cos(2.0);
    B[2][2] = m;
    count = 0;
    if (IpplInfo::noFieldCompression) {
        EXPECT_EQ(B.CompressedFraction(), 0);
        std::cout << "PASSED: Still uncompressed Field of maps" << std::endl;
    } else {
        for (lb=B.begin_if(); lb!=B.end_if(); ++lb) {
            if ( (*lb).second->IsCompressed() )
                ++count;
        }
        int reduced_count = 0;
        reduce(&count,&count+1,&reduced_count,OpAddAssign());
        EXPECT_EQ(reduced_count, 3);
    }
}

TEST(Field, FiveFields)
{
    OpalTestUtilities::SilenceTest silencer;
    int size       = 4;
    int vnodes     = 2;
    int iterations = 3;

    const unsigned Dim=3;
    Index I(size);
    Index J(size);
    Index K(size);
    FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
    Field<double,Dim> A(layout,GuardCellSizes<Dim>(2));
    Field<double,Dim> B(layout,GuardCellSizes<Dim>(2));
    Field<double,Dim> C(layout,GuardCellSizes<Dim>(2));
    Field<double,Dim> D(layout,GuardCellSizes<Dim>(2));
    Field<double,Dim> E(layout,GuardCellSizes<Dim>(1));

    FieldLayout<2> layout2(I,J,PARALLEL,PARALLEL, vnodes);

    Field<double,2> D2(layout2,GuardCellSizes<2>(2));


    A = 0.0;
    B = 0.0;
    C = 0.0;
    D = 0.0;
    E = 0.0;
    D2 = 0.0;

    A[size/2][size/2][size/2] = 512.0*(iterations + 1);
    B[size/2][size/2][size/2] = 512.0*(iterations + 2);
    C[size/2][size/2][size/2] = 512.0*(iterations + 3);
    D[size/2][size/2][size/2] = 512.0*(iterations + 4);

    for(int iter = 0 ; iter < iterations ; iter++ ) {
        std::cout << "Computing new values at iteration " << iter << " ..." << std::endl;

        E[I][J][K]  = A[I][J][K+2] + B[I][J+2][K] + C[I-2][J][K] + D[I][J][K];

        A = E;

        C = A + B;

        B = E;

        std::cout << "  iter = " << iter << ", sum(A) = " << sum(A) << std::endl;
        std::cout << "  iter = " << iter << ", sum(C) = " << sum(C) << std::endl;
        std::cout << "  iter = " << iter << ", sum(E) = " << sum(E) << std::endl;
    }

    EXPECT_NEAR(sum(A),17920,roundOffError);
    EXPECT_NEAR(sum(C),33280,roundOffError);
    EXPECT_NEAR(sum(E),17920,roundOffError);

    std::cout << A[I][J][1] << std::endl;
}

TEST(Field, Float)
{
    OpalTestUtilities::SilenceTest silencer;
    const int N=5;
    Index I(N), J(N);

    FieldLayout<2> layout2(I,J);
    Field<double,2> B(layout2),A(layout2);
    Field<float,2>  C(layout2),D(layout2);
    double d=1;
    float  f=2;
    int i=3;

    A = 1;
    B = 2.0;
    C = 3;
    D = 4.0;

    A = f;
    B = d;
    C = d;
    D = f;

    B += 2.0*A;    // d = d*d
    B += 2.0*C;    // d = d*f
    B += 2.0F*A;   // d = f*d
    B += 2.0F*C;   // d = f*f
    B += 2*A;      // d = i*d
    B += 2*C;      // d = i*f
    C += 2.0*A;    // f = d*d
    C += 2.0*C;    // f = d*f
    C += 2.0F*A;   // f = f*d
    C += 2.0F*C;   // f = f*f
    C += 2*A;      // f = f*d
    C += 2*C;      // f = f*f

    B += d*A;      // d = d*d
    B += d*C;      // d = d*f
    B += f*A;      // d = f*d
    B += f*C;      // d = f*f
    B += i*A;      // d = i*d
    B += i*C;      // d = i*f
    C += d*A;      // f = d*d
    C += d*C;      // f = d*f
    C += f*A;      // f = f*d
    C += f*C;      // f = f*f
    C += i*A;      // f = i*d
    C += i*C;      // f = i*f

    std::cout << "Results:" << std::endl;
    std::cout << "A = " << A << std::endl;
    std::cout << "B = " << B << std::endl;
    std::cout << "C = " << C << std::endl;
    EXPECT_NEAR(sum(A), N*N*2,    roundOffError);
    EXPECT_NEAR(sum(B), N*N*1129, roundOffError);
    EXPECT_NEAR(sum(C), N*N*4512, roundOffError);
}

// flyercode.cpp , Tim Williams 10/10/1997
// This is the 2D diffusion example from 1997 IPPL flyer (SC97 handout).
TEST(Field, FlyerCode)
{
    OpalTestUtilities::SilenceTest silencer;
    // Uniform 2D cartesian mesh,129x129 vertices, default spacing = (1.0,1.0):
    unsigned N = 129; Index Iverts(N); Index Jverts(N);
    typedef UniformCartesian<2> M;
    M mesh(Iverts,Jverts);

    // Boundary conditions--zero on all faces:
    Field<double,2>::bcond_container bc;
    for (int f=0; f<4; f++) bc[f] = new ConstantFace<double,2>(f,0.0);

    // Construct a 2D Field of double's representing scalar fluid field U(x,y).
    // Cell-centering on mesh, meaning 128x128 Field elements:
    CenteredFieldLayout<2,M,Cell> layout(mesh,PARALLEL,PARALLEL);
    Field<double,2,M,Cell> U(layout,GuardCellSizes<2>(1),bc);

    // Initial conditions, zero except one-cell spike in center of box:
    U = 0.0;  U[N/2][N/2] = 1000.0;
    double dt = 0.1; // Timestep

    // Compute global sum of values from field (diagnostic):
    double sumU = sum(U);
    std::cout << "sum at t = 0 is " << sumU << std::endl;
    EXPECT_NEAR(sumU, 1000, roundOffError);

    Index I(128), J(128);  // Index objects of whole-system size (all cells)

    // Run the diffusion stencil through 15 timesteps:
    for (int itime=0; itime < 15; itime++)
        U[I][J] += dt*(U[I+1][J] + U[I-1][J] - 4*U[I][J] + U[I][J+1] + U[I][J-1]);

    // Recompute global sum of values from field (diagnostic):
    std::cout << "sum at t = " << 15*dt << " is " << sum(U) << std::endl;
    EXPECT_NEAR(sumU, 1000, roundOffError);
}

TEST(Field, MinMax)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(5);
    Index J(5);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    Field<double,Dim> A(layout);
    Field<double,Dim> B(layout);
    Field<double,Dim> C(layout);

    B = 0.0;
    B[I][J] += I-2;
    C = 1.0;
    std::cout << "B" << std::endl;
    std::cout << B << std::endl;
    std::cout << sum(B) << std::endl;
    EXPECT_EQ(sum(B),0);
    std::cout << "C" << std::endl;
    std::cout << C << std::endl;
    std::cout << sum(C) << std::endl;
    EXPECT_EQ(sum(C),25);

    // Min
    A = 0.0;
    A += where( lt(B,C) , B, C );
    std::cout << "min(B,C):" << std::endl;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_EQ(sum(A),-5);

    // Max
    A = 0.0;
    A += where( gt(B,C) , B, C );
    std::cout << "max(B,C):" << std::endl;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_EQ(sum(A),30);

    // Abs
    A = 0.0;
    A += where( gt(B,0.0), B, -B);
    std::cout << "abs(B):" << std::endl;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_EQ(sum(A),30);

    // clip
    A = 0.0;
    A += where( gt(B,0.0), B, 0.0 );
    std::cout << "clip(B):" << std::endl;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_EQ(sum(A),15);
}

TEST(Field, Patches)
{
    OpalTestUtilities::SilenceTest silencer;
    const int N=5;
    Index I(N), J(N);
    BConds<double,2> bc;
    if (Ippl::getNodes() == 1) {
        bc[0] = new PeriodicFace<double,2>(0);
        bc[1] = new PeriodicFace<double,2>(1);
        bc[2] = new PeriodicFace<double,2>(2);
        bc[3] = new PeriodicFace<double,2>(3);
    }
    else {
        bc[0] = new ParallelPeriodicFace<double,2>(0);
        bc[1] = new ParallelPeriodicFace<double,2>(1);
        bc[2] = new ParallelPeriodicFace<double,2>(2);
        bc[3] = new ParallelPeriodicFace<double,2>(3);
    }

    FieldLayout<1> layout1(I);
    FieldLayout<2> layout2(I,J);
    Field<double,2> B(layout2);
    Field<double,1> C(layout1);
    Field<double,2> T2(layout2);
    Field<double,1> T1(layout1);

    int Guards = 0;
    int i,j;

    {
        Field<double,2> A(layout2,GuardCellSizes<2>(Guards),bc);

        //----------------------------------------

        assign(A[I][J] , I+J*10);
        T2 = A;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j) {
                T2[i][j] -= i+j*10.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " initializing A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[I][J] = A[J][I];
        T2 = B;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j) {
                T2[j][i] -= i+j*10.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " transposing A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[I][J] = A[I+1][J+1];
        T2 = B;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j)  {
                if ( (i<N-1)&&(j<N-1) )
                    T2[i][j] -= (i+1)+(j+1)*10.0;
                else
                    T2[i][j] += 1.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " shifting A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[4][J] = A[J][4];
        T2 = B;
        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                if ( i==4 )
                    T2[i][j] -= j+40.0;
                else
                    T2[i][j] += 1.0;
            }
        }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " copying a slice" << std::endl;

        //----------------------------------------

        C[I] = 0.1*I;
        B = -1.0;
        B[I][4] = C[I] ;
        T2 = B;
        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                if ( j==4 )
                    T2[i][j] -= i*0.1;
                else
                    T2[i][j] += 1.0;
            }
        }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " inserting a slice" << std::endl;

        //----------------------------------------

        C[I] = A[2][I] ;
        for (i=0; i<N; ++i)
            C[i] -= 2.0 + i*10.0;
        C *= C;
        EXPECT_NEAR(sum(C), 0, roundOffError);
        std::cout << " extracting a slice" << std::endl;

    }

    // Now the same tests with 1 layer of guard cells and periodic bc.
    // The answers for shifting are slightly different.
    Guards = 1;

    {
        Field<double,2> A(layout2,GuardCellSizes<2>(Guards),bc);

        //----------------------------------------

        A[I][J] = I+J*10 ;
        T2 = A;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j) {
                T2[i][j] -= i+j*10.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " initializing guarded A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[I][J] = A[J][I];
        T2 = B;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j) {
                T2[j][i] -= i+j*10.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " transposing A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[I][J] = A[I+1][J+1];
        T2 = B;
        for (i=0;i<N;++i)
            for (j=0;j<N;++j)  {
                int ii = (i+1)%N;
                int jj = (j+1)%N;
                T2[i][j] -= ii+jj*10.0;
            }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " shifting A" << std::endl;

        //----------------------------------------

        B = -1.0;
        B[4][J] = A[J][4];
        T2 = B;
        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                if ( i==4 )
                    T2[i][j] -= j+40.0;
                else
                    T2[i][j] += 1.0;
            }
        }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " copying a slice" << std::endl;

        //----------------------------------------

        C[I] = 0.1*I;
        B = -1.0;
        B[I][4] = C[I] ;
        T2 = B;
        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                if ( j==4 )
                    T2[i][j] -= i*0.1;
                else
                    T2[i][j] += 1.0;
            }
        }
        T2 *= T2;
        EXPECT_NEAR(sum(T2), 0, roundOffError);
        std::cout << " inserting a slice" << std::endl;

        //----------------------------------------

        C[I] = A[2][I] ;
        for (i=0; i<N; ++i)
            C[i] -= 2.0 + i*10.0;
        C *= C;
        EXPECT_NEAR(sum(C), 0, roundOffError);
        std::cout << " extracting a slice" << std::endl;
    }
}

TEST(Field, Reduce)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(5);
    Index J(5);
    FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL);
    Field<double,Dim> A(layout);

    Vektor<double,3> r[10];
    Vektor<double,3> rtmp[10];

    for (int i=0;i<10;i++)
        for(int d=0; d<3; d++)
            rtmp[i](d)=i*d;

    for (int i=0;i<10;i++) {
        for(int d=0; d<3; d++)
            std::cout << rtmp[i](d) << "\t";
        std::cout << std::endl;
    }

    reduce(rtmp,rtmp+10,r,OpAddAssign());

    for (int i=0;i<10;i++) {
        for(int d=0; d<3; d++)
            std::cout << r[i](d) << "\t";
        std::cout << std::endl;
    }


    A[I][J] = I+2*J+1;

    EXPECT_EQ(min(A),1);
    EXPECT_EQ(max(A),13);
    EXPECT_EQ(sum(A),175);
    EXPECT_NEAR(prod(A),3.91486e+19,1e14);
    std::cout << "min(A) = " << min(A) << std::endl;
    std::cout << "max(A) = " << max(A) << std::endl;
    std::cout << "sum(A) = " << sum(A) << std::endl;
    std::cout << "prod(A)= " << prod(A) << std::endl;

    double minv, maxv;
    minmax(A, minv, maxv);
    std::cout << "minmax(A) = (" << minv << ", " << maxv << ")" << std::endl;
    EXPECT_EQ(minv,min(A));
    EXPECT_EQ(maxv,max(A));

    std::cout << "any(A,1)= " << any(A,1.0) << std::endl;
    std::cout << "any(A,10)= " << any(A,10.0) << std::endl;
    std::cout << "any(A,1,OpGT())= " << any(A,1.0,OpGT()) << std::endl;
    std::cout << "any(A,10,OpLT())= " << any(A,10.0,OpLT()) << std::endl;
    EXPECT_EQ(any(A,1.0) ,1);
    EXPECT_EQ(any(A,10.0),1);
    EXPECT_EQ(any(A,1.0,OpGT()),1);
    EXPECT_EQ(any(A,10.0,OpLT()),1);

    NDIndex<Dim> Loc;
    double m;
    m = min(A,Loc);
    std::cout << "min(A) = " << m << " at " << Loc << std::endl;
    EXPECT_EQ(m,min(A));
    m = max(A,Loc);
    std::cout << "max(A) = " << m << " at " << Loc << std::endl;
    EXPECT_EQ(m,max(A));
}

TEST(Field, Reduceloc)
{
    Index I(5);
    Index J(5);
    FieldLayout<Dim> layout(I,J);
    Field<double,Dim> A(layout);
    Field<double,Dim> B(layout);
    Field<double,Dim> C(layout);

    A[I][J] << (I-1)*(I-1)+(J-1)*(J-1) + 1;
    NDIndex<Dim> maxloc,minloc;
    double maxval = max(A,maxloc);
    double minval = min(A,minloc);

    double known_max = 19;
    double known_min = 1;
    NDIndex<Dim> known_minloc(Index(1,1),Index(1,1));
    NDIndex<Dim> known_maxloc(Index(4,4),Index(4,4));
    EXPECT_NEAR(maxval, known_max, roundOffError);
    EXPECT_NEAR(minval, known_min, roundOffError);
    EXPECT_TRUE(maxloc == known_maxloc);
    EXPECT_TRUE(minloc == known_minloc);
}

TEST(Field, Repartition)
{
    OpalTestUtilities::SilenceTest silencer;
    const unsigned Dim=2;
    const int N=10;
    Index I(N), J(N);
    FieldLayout<Dim> layout1(I,J,PARALLEL,SERIAL);
    Field<int,Dim> A1(layout1);
    Field<int,Dim> A2(layout1);

    assign(A1[I][J], I + 10*J);
    assign(A2[I][J], J + 10*I);

    // an alternative way to generate the list of local vnode domains,
    // directly from a FieldLayout
    std::vector< NDIndex<Dim> > newDomains2;

    FieldLayout<Dim>       layout2(I,J,SERIAL,PARALLEL);

    FieldLayout<Dim>::iterator_iv lvnode, endvnode=layout2.end_iv();
    for (lvnode = layout2.begin_iv(); lvnode != endvnode; ++lvnode)
        newDomains2.push_back(NDIndex<Dim>((*(*lvnode).second).getDomain()));

    layout1.Repartition(&newDomains2[0],&newDomains2[0]+newDomains2.size());

    std::cout << "My portion of Field A1 = " << A1 << std::endl;
    std::cout << "My portion of Field A2 = " << A2 << std::endl;
    EXPECT_NEAR(sum(A1),4950,roundOffError);
    EXPECT_NEAR(sum(A2),4950,roundOffError);
}

TEST(Field, ScalarIndexing)
{
    const unsigned nx = 4, ny = 4, nz = 4;
    Index I(nx), J(ny), K(nz);
    FieldLayout<D3> layout(I,J,K);

    // Instantiate and initialize scalar, vector, tensor fields:
    Field<double,D3> scalarFld(layout);
    double scalar = 1.0;
    scalarFld << scalar;
    Field<Vektor<double,D3>,D3> vectorFld(layout);
    Vektor<double, D3> vector(1.0,2.0,3.0);
    vectorFld << vector;
    Field<Tenzor<double,D3>,D3,UniformCartesian<D3> > tensorFld(layout);
    Tenzor<double, D3> tensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    tensorFld << tensor;
    Field<SymTenzor<double,D3>,D3,UniformCartesian<D3> > symTensorFld(layout);
    SymTenzor<double, D3> symTensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    symTensorFld << symTensor;

    // Now try the scalar indexing:
    double scalar1 = 0.0;
    Vektor<double, D3> vector1(0.0, 0.0, 0.0);
    Tenzor<double, D3> tensor1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    SymTenzor<double, D3> symTensor1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    scalar1 = scalarFld[1][1][1].get();
    vector1 = vectorFld[1][1][1].get();
    tensor1 = tensorFld[1][1][1].get();
    symTensor1 = symTensorFld[1][1][1].get();

    EXPECT_TRUE((scalar1 == scalar));
    EXPECT_TRUE((vector1 == vector));
    EXPECT_TRUE((tensor1 == tensor));
    EXPECT_TRUE((symTensor1 == symTensor));
}

TEST(Field, Single)
{
    OpalTestUtilities::SilenceTest silencer;
    const int Dim=2;
    int n=5;
    Index I(n),J(n);
    FieldLayout<Dim> layout(I,J);
    Field<double,Dim> A(layout);
    int i,j;

    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
            A[i][j] = i+10*j;
        }
    }
    for (j=n-1; j>=0; --j) {
        for (i=0; i<n; ++i) {
            EXPECT_NEAR(A[i][j].get(),i+10*j,roundOffError);
        }
    }
}

TEST(Field, SimpleTest1)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(5);
    Index J(5);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    Field<double,Dim> A(layout);
    Field<int,Dim> B(layout,GuardCellSizes<Dim>(1));
    Field<int,Dim>::iterator p;
    int i;

    A = -1.0 ;
    for (p=B.begin(), i=0; p!=B.end(); ++p, ++i) *p = i;

    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),-25,roundOffError);
    std::cout << B << std::endl;
    std::cout << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),600,roundOffError);

    A[I][J] = B[I-1][J+1] + B[I+1][J-1] ;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),1200,roundOffError);

    B.fillGuardCells();
    A[I][J] = B[I-1][J+1] + B[I+1][J-1];
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),1200,roundOffError);
}

TEST(Field, SimpleTest2)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(4);
    Index J(4);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));
    Field<double,Dim>::iterator p;
    int i;

    A.fillGuardCells();

    for (p=A.begin(), i=0; p!=A.end(); ++p, ++i) *p = i;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),240,roundOffError);
}

TEST(Field, SimpleTest3)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(4);
    Index J(4);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    Field<double,Dim> A(layout);
    Field<double,Dim> B(layout,GuardCellSizes<Dim>(1));
    Field<double,Dim>::iterator p;
    int i;

    A = -1.0 ;
    for (p=B.begin(), i=0; p!=B.end(); ++p, ++i) *p = i+1;
    B = B*2.0;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),-16,roundOffError);
    std::cout << B << std::endl;
    std::cout << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),512,roundOffError);

    A[I][J] = (B[I+1][J+1]+B[I+1][J-1]+B[I-1][J+1]+B[I-1][J-1])/8.0;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),256,roundOffError);
    // t4
    A[I][J]= (I-2.5)*(I-1.5) + (J-1.5)*(J-2.5) ;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),40,roundOffError);
}

TEST(Field, SimpleTest5)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(4);
    Index J(4);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);

    Field<double,Dim> A(layout);
    Field<double,Dim> B(layout);

    A = 0;

    A[I][J] += (I-2)*(I-2) + (J-2)*(J-2);
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),48,roundOffError);

    A -= 1.0;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),32,roundOffError);

    B[I][J] = I;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),32,roundOffError);
    std::cout << B << std::endl;
    std::cout << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),24,roundOffError);

    A *= B;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),28,roundOffError);

    //t6
    A = 0;
    Index I1(1,3),J1(1,3);
    A[I1][J1] += I1 + J1/100.0;
    std::cout << A << std::endl;
    std::cout << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),18.18,roundOffError);
}

TEST(Field, SimpleTest7)
{
    OpalTestUtilities::SilenceTest silencer;
    Index I(10);       // Index on [0..9]
    Index Low(5);      // Index on [0..4]
    Index High(5,9);   // Index on [5..9]
    Index Even(0,9,2); // Index on [0..9] stride 2
    Index Odd(1,9,2);  // Index on [1..9] stride 2
    Index Edge(0,9,9); // Index on [0..9] stride 9; the edge cells

    NDIndex<D1> domain;
    domain[0] = I;
    FieldLayout<D1> layout(domain);
    Field<double,D1> A(layout);

    A = 0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),0,roundOffError);

    A[I-5] = 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),5,roundOffError);

    A[I+5] = 2.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),15,roundOffError);

    A[I*2] -= 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),10,roundOffError);

    A[I*2+1] -= 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),5,roundOffError);

    A[I*9] = 3.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),10,roundOffError);

    std::cout << "Should be the same as :" << std::endl;

    A = 0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),0,roundOffError);

    A[Low] = 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),5,roundOffError);

    A[High] = 2.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),15,roundOffError);

    A[Even] -= 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),10,roundOffError);

    A[Odd] -= 1.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),5,roundOffError);

    A[Edge] = 3.0;

    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),10,roundOffError);
}

TEST(Field, SimpleTest8)
{
    OpalTestUtilities::SilenceTest silencer;
    Vektor<double,D3> boxMin(2,3,4);
    Vektor<double,D3> boxMax(5,6,7);
    Vektor<double,D3> h(1.0e-2, 1.0e-2, 0.5e-2);

    Vektor<double,D3> p1(8,9,10);

    Vektor<unsigned,D3> N(static_cast<unsigned> (ceil( (abs(boxMin[0])+boxMax[0])/h[0])),
                           static_cast<unsigned> (ceil( (abs(boxMin[1])+boxMax[1])/h[1])),
                           static_cast<unsigned> (ceil( (abs(boxMin[2])+boxMax[2])/h[2])));

    std::cout << "orig= " << boxMin << " maxext= " << boxMax << std::endl;
    std::cout << "h=  " << h << " N= " << N << std::endl;
    std::cout << "p1= " << p1 << std::endl;

    /*
      p1 -> n
    */
    Vektor<double, D3> np1 = -(boxMin-p1)/h;
    std::cout << "n(p1)= " << np1 << std::endl;
    EXPECT_NEAR(np1[0],600,roundOffError);
    EXPECT_NEAR(np1[1],600,roundOffError);
    EXPECT_NEAR(np1[2],1200,roundOffError);
}

TEST(Field, Subdivide)
{
    OpalTestUtilities::SilenceTest silencer;
    const int Dim=2;
    int N=5, i, j;
    Index II(N+2),JJ(N+2);
    Index I(0,N),J(0,N);
    FieldLayout<Dim> layout(II,JJ,PARALLEL,PARALLEL,4);
    Field<double,Dim> A(layout),B(layout);

    A=-1.0;
    assign(B[II][JJ], II+JJ*10);
    A[I][J] = B[J+1][I+1];

    for (j=0; j<=N+1; ++j) {
        for(i=0; i<=N+1; ++i)
            std::cout << " " << B[i][j];
        std::cout << std::endl;
    }
    std::cout << std::endl;
    EXPECT_NEAR(sum(B),1617,roundOffError);

    for (Field<double,Dim>::iterator p=B.begin(); p!=B.end(); ++p)
        std::cout << " " << *p;
    std::cout << std::endl;

    for (j=0; j<=N+1; ++j) {
        for(i=0; i<=N+1; ++i)
            std::cout << " " << A[i][j].get();
        std::cout << std::endl;
    }
    std::cout << std::endl;
    EXPECT_NEAR(sum(A),1373,roundOffError);
}

TEST(Field, SubReadTest)
{
    OpalTestUtilities::SilenceTest silencer;
    int sizeX   = 3;
    int sizeY   = 3;
    int vnodes  = 2;
    int vnodes2 = 2;
    int fields  = 4;

    const unsigned Dim=2;
    typedef int T;
    Index I(1, sizeX), I2(1 + sizeX/3, 1 + 2*sizeX/3), I3(0,2*sizeX + 1);
    Index J(1, sizeY), J2(1 + sizeY/3, 1 + 2*sizeY/3), J3(0,2*sizeY + 1);
    NDIndex<Dim> domain(I, J);
    NDIndex<Dim> subdom(I2, J2);
    NDIndex<Dim> bigdom(I3, J3);
    FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL, vnodes);
    FieldLayout<Dim> layout2(I,J,PARALLEL,PARALLEL, vnodes2);
    FieldLayout<Dim> layout3(I3,J3,PARALLEL,PARALLEL, vnodes2);
    Field<T,Dim> A(layout);
    Field<T,Dim> A2(layout);
    Field<T,Dim> B(layout2);
    Field<T,Dim> B2(layout);
    Field<T,Dim> C(layout3);
    Field<T,Dim> C2(layout3);

    // initialize data
    std::cout << "Initializing A ..." << std::endl;
    A[I][J] = I + 10*(J-1);
    //FieldPrint<T,Dim> fp(A);
    //FieldPrint<T,Dim> fpb(B);
    //std::cout << "Initial A:" << std::endl;
    //fp.print(domain);
    //std::cout << "A subset:" << std::endl;
    //fp.print(subdom);

    // write A to disk N times
    std::cout << "====================== writing ====================" << std::endl;
    DiscField<Dim> dcf("ackfiledata", "ackfiledata.config", fields, "data");
    std::cout << "Writing field " << 0 << " ..." << std::endl;
    dcf.write(A, 0);
    A2 = A;
    for (int i=1; i < fields; ++i) {
        std::cout << "Incrementing field ..." << std::endl;
        A2 += 100;

        std::cout << "Writing field " << i << " ..." << std::endl;
        dcf.write(A2, i);
    }

    // read in the first and last field
    std::cout << "====================== reading ====================" << std::endl;
    DiscField<Dim> dcf2("ackfiledata", "ackfiledata.config");
    std::cout << "Reading in field 0 ..." << std::endl;
    dcf2.read(B, 0);
    //std::cout << "Field just read: " << std::endl;
    //fpb.print(domain);
    B2 = B;
    B2 -= A;
    NDIndex<Dim> mloc;
    std::cout << "Checking field just read: min and max of diff should be zero." << std::endl;
    std::cout << "     min(B-A) = " << min(B2, mloc) << std::endl;
    std::cout << "  minloc(B-A) = " << mloc << std::endl;
    std::cout << "     max(B-A) = " << max(B2, mloc) << std::endl;
    std::cout << "  maxloc(B-A) = " << mloc << std::endl;
    std::cout << "       sum(A) = " << sum(A) << std::endl;
    std::cout << "       sum(B) = " << sum(B) << std::endl;
    EXPECT_NEAR(min(B2, mloc),0,roundOffError);
    EXPECT_NEAR(max(B2, mloc),0,roundOffError);
    EXPECT_NEAR(sum(A),108,roundOffError);
    EXPECT_NEAR(sum(B),108,roundOffError);
    A2 = 0;
    A2[I2][J2] = I2 + 10*(J2-1);
    std::cout << "Reading in field 0 again, with subdomain " << subdom <<" ..."<<std::endl;
    B = 0;
    dcf2.read(B, subdom, 0);
    //std::cout << "Field just read: " << std::endl;
    //fpb.print(domain);
    B2 = B;
    B2 -= A2;
    std::cout << "Checking field just read: min and max of diff should be zero." << std::endl;
    std::cout << "     min(B-A2) = " << min(B2, mloc) << std::endl;
    std::cout << "  minloc(B-A2) = " << mloc << std::endl;
    std::cout << "     max(B-A2) = " << max(B2, mloc) << std::endl;
    std::cout << "  maxloc(B-A2) = " << mloc << std::endl;
    std::cout << "       sum(A2) = " << sum(A2) << std::endl;
    std::cout << "        sum(B) = " << sum(B) << std::endl;
    EXPECT_NEAR(min(B2, mloc),0,roundOffError);
    EXPECT_NEAR(max(B2, mloc),0,roundOffError);
    EXPECT_NEAR(sum(A2),70,roundOffError);
    EXPECT_NEAR(sum(B), 70,roundOffError);

    std::cout << "Reading in field 0 again into big array, with subdomain " << subdom <<" ..."<<std::endl;
    C = 0;
    C2 = 0;
    C2[I2][J2] = I2 + 10*(J2-1);
    dcf2.read(C, subdom, 0);
    std::cout << "Checking field just read: min and max of diff should be zero." << std::endl;
    std::cout << "     min(C2-C) = " << min(C2-C, mloc) << std::endl;
    std::cout << "  minloc(C2-C) = " << mloc << std::endl;
    std::cout << "     max(C2-C) = " << max(C2-C, mloc) << std::endl;
    std::cout << "  maxloc(C2-C) = " << mloc << std::endl;
    std::cout << "       sum(C2) = " << sum(C2) << std::endl;
    std::cout << "        sum(C) = " << sum(C) << std::endl;
    EXPECT_NEAR(min(C2-C, mloc),0,roundOffError);
    EXPECT_NEAR(max(C2-C, mloc),0,roundOffError);
    EXPECT_NEAR(sum(C2),70,roundOffError);
    EXPECT_NEAR(sum(C), 70,roundOffError);
}

TEST(Field, Transpose)
{
    OpalTestUtilities::SilenceTest silencer;
    const int N=10;

    Index I(N),J(N),K(N),L(N);
    NDIndex<D2> IJ(I,J);
    NDIndex<D3> IJK(IJ,K);
    NDIndex<D4> IJKL(IJK,L);
    e_dim_tag dist[D4] = { PARALLEL,PARALLEL,SERIAL,SERIAL };
    FieldLayout<D2> layout2(IJ,dist);
    FieldLayout<D3> layout3(IJK,dist);
    FieldLayout<D4> layout4(IJKL,dist);

    Field<int,D2> A1(layout2),A2(layout2);
    Field<int,D3> B1(layout3),B2(layout3);
    Field<int,D4> C1(layout4),C2(layout4);
    int ii = 1;
    int jj = 10;
    int kk = 100;
    int ll = 1000;
    int s;

    A1[I][J] = ii*I + jj*J;
    A2[J][I] = A1[I][J] ;
    A2[I][J] -= ii*J + jj*I ;
    A2 *= A2;
    s = sum(A2);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "2D General transpose: " << (s?"FAILED":"PASSED") << std::endl;

    B1[I][J][K] = ii*I + jj*J + kk*K;
    B2[I][J][K] = B1[I][K][J] ;
    B2[I][J][K] -= ii*I + jj*K + kk*J;
    B2 *= B2;
    s = sum(B2);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "3D General transpose: " << (s?"FAILED":"PASSED") << std::endl;

    C1[I][J][K][L] = ii*I + jj*J + kk*K + ll*L ;
    C2[I][J][L][K] = C1[I][J][K][L] ;
    C2[I][J][K][L] -= ii*I + jj*J + kk*L + ll*K ;
    C2 *= C2;
    s = sum(C2);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "4D serial transpose: "<< (s?"FAILED":"PASSED") << std::endl;

    int a,b;
    s = 0;
    for (a=0; a<N; ++a)
        {
            A1[I][J] = B1[a][I][J];
            A1[I][J] -= a*ii + I*jj + J*kk;
            A1 *= A1;
            s += sum(A1);
        }
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << std::endl;

    s = 0;
    for (a=0; a<N; ++a)
        {
            A1[I][J] = B1[I][a][J];
            A1[I][J] -= I*ii + a*jj + J*kk;
            A1 *= A1;
            s += sum(A1);
        }
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << std::endl;

    s = 0;
    for (a=0; a<N; ++a)
        {
            A1[I][J] = B1[I][J][a];
            A1[I][J] -= I*ii + J*jj + a*kk;
            A1 *= A1;
            s += sum(A1);
        }
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "Serial 2 from 3 " << (s ? "FAILED " : "PASSED ") << std::endl;

    A1[I][J] = ii*I + jj*J;
    for (a=0; a<N; ++a)
        B1[a][I][J] = A1[I][J] ;
    B1[K][I][J] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << std::endl;

    for (a=0; a<N; ++a)
        B1[I][a][J] = A1[I][J] ;
    B1[I][K][J] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << std::endl;

    for (a=0; a<N; ++a)
        B1[I][J][a] = A1[I][J] ;
    B1[I][J][K] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "Serial 3 from 2 " << (s ? "FAILED " : "PASSED ") << std::endl;

    for (a=0; a<N; ++a)
        for (b=0; b<N; ++b)
            C1[I][J][a][b] = A1[I][J] ;
    C1[I][J][K][L] -= I*ii + J*jj;
    C1 *= C1;
    s = sum(C1);
    EXPECT_NEAR(s, 0, roundOffError);
    std::cout << "Serial 4 from 2 " << (s ? "FAILED " : "PASSED ") << std::endl;
}

namespace {
    void check( Field<int,D3>& f, int s1, int s2, int s3, int /*test*/)
    {
        Index I = f.getIndex(0);
        Index J = f.getIndex(1);
        Index K = f.getIndex(2);
        f[I][J][K] -= s1*I + s2*J + s3*K;
        int sum_f = sum(f);
        EXPECT_NEAR(sum_f, 0, roundOffError);
    }
}

TEST(Field, Transpose2)
{
    const int N=10;

    Index I(N),J(N),K(N);
    FieldLayout<D3> layout_ppp(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);
    FieldLayout<D3> layout_spp(I,J,K,SERIAL,PARALLEL,PARALLEL,8);
    FieldLayout<D3> layout_psp(I,J,K,PARALLEL,SERIAL,PARALLEL,8);
    FieldLayout<D3> layout_pps(I,J,K,PARALLEL,PARALLEL,SERIAL,8);

    Field<int,D3> A(layout_ppp);
    Field<int,D3> B(layout_spp);
    Field<int,D3> C(layout_psp);
    Field<int,D3> D(layout_pps);
    int ii = 1;
    int jj = 10;
    int kk = 100;

    assign( A[I][J][K] , ii*I + jj*J + kk*K );
    B = A;
    C = A;
    D = A;
    check(B,ii,jj,kk,1);
    check(C,ii,jj,kk,2);
    check(D,ii,jj,kk,3);

    B[I][J][K] = A[I][J][K];
    C[I][J][K] = A[I][J][K];
    D[I][J][K] = A[I][J][K];
    check(B,ii,jj,kk,4);
    check(C,ii,jj,kk,5);
    check(D,ii,jj,kk,6);

    B[I][K][J] = A[I][J][K];
    C[I][K][J] = A[I][J][K];
    D[I][K][J] = A[I][J][K];
    check(B,ii,kk,jj,7);
    check(C,ii,kk,jj,8);
    check(D,ii,kk,jj,9);

    B[J][I][K] = A[I][J][K];
    C[J][I][K] = A[I][J][K];
    D[J][I][K] = A[I][J][K];
    check(B,jj,ii,kk,10);
    check(C,jj,ii,kk,11);
    check(D,jj,ii,kk,12);

    B[J][K][I] = A[I][J][K];
    C[J][K][I] = A[I][J][K];
    D[J][K][I] = A[I][J][K];
    check(B,jj,kk,ii,13);
    check(C,jj,kk,ii,14);
    check(D,jj,kk,ii,15);

    B[K][I][J] = A[I][J][K];
    C[K][I][J] = A[I][J][K];
    D[K][I][J] = A[I][J][K];
    check(B,kk,ii,jj,16);
    check(C,kk,ii,jj,17);
    check(D,kk,ii,jj,18);

    B[K][J][I] = A[I][J][K];
    C[K][J][I] = A[I][J][K];
    D[K][J][I] = A[I][J][K];
    check(B,kk,jj,ii,19);
    check(C,kk,jj,ii,20);
    check(D,kk,jj,ii,21);
}