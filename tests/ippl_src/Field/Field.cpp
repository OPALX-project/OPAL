#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Field/BareField.h"
#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Index/Index.h"
#include "Meshes/UniformCartesian.h"

#include <map>
#include <iostream>

constexpr unsigned Dim = 2;
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
            // testmsg << "PASSED: Uncompressed one after assigning in guard cell"
        }
        else {
            EXPECT_EQ(reduced_count,2);
            //testmsg << "PASSED: Uncompressed two after assigning in guard cell"
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
        //testmsg << "FAILED: Compressed somewhere after squaring" << endl;
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
            //testmsg << "PASSED: Uncompressed one after squaring"
        }
        else {
            EXPECT_EQ(reduced_count,2);
            //testmsg << "PASSED: Uncompressed two after squaring"
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
        //  testmsg << "PASSED: Still uncompressed Field of maps" << endl;
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

TEST(Field, Patches)
{
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
        // testmsg << " initializing A" << endl;

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
        //testmsg << " transposing A" << endl;

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
        // testmsg << " shifting A" << endl;

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
        //testmsg << " copying a slice" << endl;

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
        //testmsg << " inserting a slice" << endl;

        //----------------------------------------

        C[I] = A[2][I] ;
        for (i=0; i<N; ++i)
            C[i] -= 2.0 + i*10.0;
        C *= C;
        EXPECT_NEAR(sum(C), 0, roundOffError);
        //testmsg << " extracting a slice" << endl;

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
        //testmsg << " initializing guarded A" << endl;

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
        //testmsg << " transposing A" << endl;

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
        //testmsg << " shifting A" << endl;

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
        //testmsg << " copying a slice" << endl;

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
        //testmsg << " inserting a slice" << endl;

        //----------------------------------------

        C[I] = A[2][I] ;
        for (i=0; i<N; ++i)
            C[i] -= 2.0 + i*10.0;
        C *= C;
        EXPECT_NEAR(sum(C), 0, roundOffError);
        //testmsg << " extracting a slice" << endl;
    }
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

TEST(Field, Transpose)
{
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
    //testmsg << "2D General transpose: " << (s?"FAILED":"PASSED") << endl;

    B1[I][J][K] = ii*I + jj*J + kk*K;
    B2[I][J][K] = B1[I][K][J] ;
    B2[I][J][K] -= ii*I + jj*K + kk*J;
    B2 *= B2;
    s = sum(B2);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "3D General transpose: " << (s?"FAILED":"PASSED") << endl;

    C1[I][J][K][L] = ii*I + jj*J + kk*K + ll*L ;
    C2[I][J][L][K] = C1[I][J][K][L] ;
    C2[I][J][K][L] -= ii*I + jj*J + kk*L + ll*K ;
    C2 *= C2;
    s = sum(C2);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "4D serial transpose: "<< (s?"FAILED":"PASSED") << endl;

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
    //testmsg << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << endl;

    s = 0;
    for (a=0; a<N; ++a)
        {
            A1[I][J] = B1[I][a][J];
            A1[I][J] -= I*ii + a*jj + J*kk;
            A1 *= A1;
            s += sum(A1);
        }
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << endl;

    s = 0;
    for (a=0; a<N; ++a)
        {
            A1[I][J] = B1[I][J][a];
            A1[I][J] -= I*ii + J*jj + a*kk;
            A1 *= A1;
            s += sum(A1);
        }
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "Serial 2 from 3 " << (s ? "FAILED " : "PASSED ") << endl;

    A1[I][J] = ii*I + jj*J;
    for (a=0; a<N; ++a)
        B1[a][I][J] = A1[I][J] ;
    B1[K][I][J] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

    for (a=0; a<N; ++a)
        B1[I][a][J] = A1[I][J] ;
    B1[I][K][J] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

    for (a=0; a<N; ++a)
        B1[I][J][a] = A1[I][J] ;
    B1[I][J][K] -= I*ii + J*jj;
    B1 *= B1;
    s = sum(B1);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "Serial 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

    for (a=0; a<N; ++a)
        for (b=0; b<N; ++b)
            C1[I][J][a][b] = A1[I][J] ;
    C1[I][J][K][L] -= I*ii + J*jj;
    C1 *= C1;
    s = sum(C1);
    EXPECT_NEAR(s, 0, roundOffError);
    //testmsg << "Serial 4 from 2 " << (s ? "FAILED " : "PASSED ") << endl;
}

namespace {
    void check( Field<int,D3>& f, int s1, int s2, int s3, int test)
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