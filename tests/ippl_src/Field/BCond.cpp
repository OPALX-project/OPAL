#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "AppTypes/Vektor.h"
#include "Field/Field.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "FieldLayout/FieldLayout.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/CartesianCentering.h"
#include "Utility/FieldDebug.h"

#include <fstream>
#include <iostream>
#include <string>

namespace {
    void hardCodedOutput(std::string filename); // Prototype of function defined below.
    bool thediff(std::string filename1, std::string filename2);
}

constexpr unsigned Dim = 2;
constexpr double   roundOffError = 1e-10;

// template definition
CenteringEnum CCCEnums<2U,2U,0U>::vectorFace[2U*2U];

TEST(Field, Left)
{
    OpalTestUtilities::SilenceTest silencer;

    bool b1=true; bool b2=false;
    std::cout << "bool output key: true = " << b1 << " ; false = " << b2 << std::endl;

    // Sizes
    const unsigned D = 1;
    int ncells = 7, nverts = 8;
    int vnodes = 2;

    // Mesh
    NDIndex<1U> verts;
    for (unsigned int d=0; d<D; d++) verts[d] = Index(nverts);
    typedef UniformCartesian<1U> M1;
    M1 mesh(verts);

    // Layouts
    CenteredFieldLayout<1U,M1,Cell> layoutCell(mesh,PARALLEL,vnodes);
    CenteredFieldLayout<1U,M1,Vert> layoutVert(mesh,PARALLEL,vnodes);

    // Guard cells
    GuardCellSizes<1U> gc(1);

    // Boundary conditions
    BConds<double,1U,M1,Cell> cbc;
    BConds<double,1U,M1,Vert> vbc;
    for (unsigned int face=0; face<2*1U; face++) {
        cbc[face] = new NegReflectFace<double,1U,M1,Cell>(face);
        vbc[face] = new NegReflectFace<double,1U,M1,Vert>(face);
    }

    // Fields
    Field<double,1U,M1,Cell> A(layoutCell,gc,cbc);
    Field<double,1U,M1,Vert> B(layoutVert,gc,vbc);
    Field<double,1U,M1,Cell> A0(layoutCell,gc,cbc);
    Field<double,1U,M1,Vert> B0(layoutVert,gc,vbc);

    // Initial values (duplicate in A0,B0 for FieldDebug output w/o changing
    // dirty_m of actual A and B Fields):
    Index I(ncells);
    Index Iv(nverts);
    A[I] = I;
    A0[I] = I;
    B[Iv] = Iv;
    B0[Iv] = Iv;

    // Output State of dirty flags prior to "stencil" assignment:
    std::cout << "!!!!!!!!!!!!! BEFORE !!!!!!!!!!!!!" << std::endl;
    std::cout << "B.isDirty() = " << B.isDirty() << " ; "
         << "A.isDirty() = " << A.isDirty() << std::endl;
    EXPECT_FALSE(A.isDirty());
    EXPECT_FALSE(B.isDirty());

    // Output (copies of) initial values):
    std::cout << "[[[[[[[ A ]]]]]]]" << std::endl;
    fp1(A0);
    EXPECT_NEAR(sum(A),21,roundOffError);
    std::cout << "[[[[[[[ B ]]]]]]]" << std::endl;
    fp1(B0);
    EXPECT_NEAR(sum(B),28,roundOffError);

    // "Stencil" assignment:
    B[I + 1] = B[I + 1] + A[I + 1];

    // Output State of dirty flags after to "stencil" assignment:
    std::cout << "!!!!!!!!!!!!! AFTER !!!!!!!!!!!!!" << std::endl;
    std::cout << "B.isDirty() = " << B.isDirty() << " ; "
         << "A.isDirty() = " << A.isDirty() << std::endl;
    EXPECT_FALSE(A.isDirty());
    EXPECT_FALSE(B.isDirty());

    // Output resulting value of B:
    std::cout << "[[[[[[[ B ]]]]]]]" << std::endl;
    fp1(B);
    EXPECT_NEAR(sum(B),43,roundOffError);
}

TEST(Field, BCSimple)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(5);
    Index J(5);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    typedef UniformCartesian<Dim> M;
    typedef Cell C;
    Field<double,Dim,M,C> B(layout);

    // Set initial boundary conditions.
    BConds<double,Dim,M,C> bc;
    if (Ippl::getNodes() == 1) {
        bc[0] = new PeriodicFace<double,Dim,M,C>(2);
        bc[1] = new PeriodicFace<double,Dim,M,C>(3);
    }
    else {
        bc[0] = new ParallelPeriodicFace<double,Dim,M,C>(2);
        bc[1] = new ParallelPeriodicFace<double,Dim,M,C>(3);
    }
    bc[2] = new PosReflectFace<double,Dim,M,C>(0);
    bc[3] = new ZeroFace<double,Dim,M,C>(1);

    // An array for testing.
    Field<double,Dim,M,C> A(layout,GuardCellSizes<Dim>(2),bc);

    // Override one.
    A.getBConds()[2] = new NegReflectFace<double,Dim,M,C>(0);

    A = 0.0;
    A[I][J] += I*10 + J;
    std::cout << "A = " << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A),550,roundOffError);
    B[I][J] = A[I-2][J-2];
    std::cout << "B = " << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),110,roundOffError);
    B[I][J] = A[I+2][J+2];
    std::cout << "B = " << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),480,roundOffError);
    B[I][J] = A[I+2][J-2];
    std::cout << "B = " << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),480,roundOffError);
    B[I][J] = A[I-2][J+2];
    std::cout << "B = " << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B),110,roundOffError);
}

namespace {
    // define some helper functions for computing projections
    Vektor<double,Dim> proj1(const Vektor<double,Dim> &v) {
        Vektor<double,Dim> vtmp = v;
        vtmp[1] = -vtmp[1];
        return vtmp;
    }

    Vektor<double,Dim> proj2(const Vektor<double,Dim> &v) {
        Vektor<double,Dim> vtmp = v;
        vtmp[0] = -vtmp[0];
        return vtmp;
    }

    double componentProj1(double vc) {
        vc = - vc;
        return vc;
    }

    double componentProj2(double vc) {
        vc = - vc;
        return vc;
    }
}

TEST(Field, BCSimple2)
{
    OpalTestUtilities::SilenceTest silencer;

    Index I(5),J(5);
    FieldLayout<Dim> layout(I,J);
    typedef UniformCartesian<Dim> M;
    typedef Cell C;
    Field< Vektor<double,Dim>, Dim, M, C > B(layout);

    // Set initial boundary conditions.
    BConds< Vektor<double,Dim>, Dim, M, C > bc;
    if (Ippl::getNodes() == 1) {
        bc[0] = new PeriodicFace< Vektor<double,Dim>, Dim, M, C >(0);
        bc[1] = new PeriodicFace< Vektor<double,Dim>, Dim, M, C >(1);
    }
    else {
        bc[0] = new ParallelPeriodicFace< Vektor<double,Dim>, Dim, M, C >(0);
        bc[1] = new ParallelPeriodicFace< Vektor<double,Dim>, Dim, M, C >(1);
    }
    bc[2] = new FunctionFace< Vektor<double,Dim>, Dim, M, C >(proj1,2);
    bc[3] = new FunctionFace< Vektor<double,Dim>, Dim, M, C >(proj2,3);

    // construct a FieldSpec object
    FieldSpec< Vektor<double,Dim>, Dim, M, C >
        Spec(layout,GuardCellSizes<Dim>(1),bc);
    FieldSpec< Vektor<double,Dim>, Dim, M, C > Sp1(Spec);
    FieldSpec< Vektor<double,Dim>, Dim, M, C > Sp2(layout);
    Sp2.set_BC(bc);
    Sp2.set_GC(GuardCellSizes<Dim>(1));
    Sp2 = Sp1;

    // An array for testing.
    Field< Vektor<double,Dim>, Dim, M, C > A(Spec);

    Field< Vektor<double,Dim>, Dim, M, C>::iterator ip;

    int i = 0, j = 0;

    A.Uncompress(); // Needed because of bug in using iterator (TJW 3/31/1997)
    for( ip = A.begin() ; ip != A.end() ; ++ip ) {
        *ip = i + j*10;
        i++;
        if( (i % 5) == 0 ) {
            i = 0;
            j++;
        }
    }
    A.fillGuardCells();
    std::cout << A << std::endl;
    std::cout << "sum A = " << sum(A) << std::endl;
    EXPECT_NEAR(sum(A)[0],550,roundOffError);
    EXPECT_NEAR(sum(A)[1],550,roundOffError);
    // A.writeb("atest");
    assign(B[I][J], A[I-1][J-1]);
    std::cout << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B)[0],350,roundOffError);
    EXPECT_NEAR(sum(B)[1],330,roundOffError);
    // B.writeb("btest1");
    assign(B[I][J], A[I+1][J+1]);
    std::cout << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B)[0],330,roundOffError);
    EXPECT_NEAR(sum(B)[1],750,roundOffError);
    // B.writeb("btest2");

    // Now try same using componentwise specification of y BC:
    // Set initial boundary conditions.
    BConds< Vektor<double,Dim>, Dim, M, C > bc2;
    if (Ippl::getNodes() == 1) {
        bc2[0] = new PeriodicFace< Vektor<double,Dim>, Dim, M, C >(0);
        bc2[1] = new PeriodicFace< Vektor<double,Dim>, Dim, M, C >(1);
    }
    else {
        bc2[0] = new ParallelPeriodicFace< Vektor<double,Dim>, Dim, M, C >(0);
        bc2[1] = new ParallelPeriodicFace< Vektor<double,Dim>, Dim, M, C >(1);
    }
    bc2[2] = new ComponentFunctionFace< Vektor<double,Dim>, Dim, M, C >
        (componentProj1,2,0);
    bc2[3] = new ComponentFunctionFace< Vektor<double,Dim>, Dim, M, C >
        (componentProj2,3,1);
    bc2[2] = new ComponentFunctionFace< Vektor<double,Dim>, Dim, M, C >
        (componentProj1,2,0);
    bc2[3] = new ComponentFunctionFace< Vektor<double,Dim>, Dim, M, C >
        (componentProj2,3,1);
    // Construct a FieldSpec object:
    FieldSpec< Vektor<double,Dim>, Dim, M, C >
        Spec2(layout,GuardCellSizes<Dim>(1),bc2);
    // Another Field for testing:
    Field< Vektor<double,Dim>, Dim, M, C > A2(Spec2);
    i = 0; j = 0;
    A2.Uncompress(); // Needed because of bug in using iterator (TJW 3/31/1997)
    for( ip = A2.begin() ; ip != A2.end() ; ++ip ) {
        *ip = i + j*10;
        i++;
        if( (i % 5) == 0 ) {
            i = 0;
            j++;
        }
    }
    A2.fillGuardCells();
    std::cout << A2 << std::endl;
    std::cout << "sum A2 = " << sum(A2) << std::endl;
    EXPECT_NEAR(sum(A2)[0],550,roundOffError);
    EXPECT_NEAR(sum(A2)[1],550,roundOffError);
    // A.writeb("a2test");
    assign(B[I][J], A[I-1][J-1]);
    std::cout << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B)[0],350,roundOffError);
    EXPECT_NEAR(sum(B)[1],330,roundOffError);
    // B.writeb("b2test1");
    assign(B[I][J], A[I+1][J+1]);
    std::cout << B << std::endl;
    std::cout << "sum B = " << sum(B) << std::endl;
    EXPECT_NEAR(sum(B)[0],330,roundOffError);
    EXPECT_NEAR(sum(B)[1],750,roundOffError);
    // B.writeb("b2test2");
}

TEST(Field, BC)
{
    // For writing file output to compare against hardcoded correct file output:
    Inform fdi(NULL,"text.test.TestBC",Inform::OVERWRITE,0);
    setInform(fdi);

    Index I(5);
    Index J(5);
    NDIndex<Dim> domain;
    domain[0] = I;
    domain[1] = J;
    FieldLayout<Dim> layout(domain);
    typedef UniformCartesian<Dim> M;

    // Set Cell-centered boundary conditions.
    BConds<double,Dim,M,Cell> cbc;
    cbc[0] = new NegReflectFace<double,Dim,M,Cell>(0);
    cbc[1] = new ZeroFace<double,Dim,M,Cell>(1);
    cbc[2] = new ParallelPeriodicFace<double,Dim,M,Cell>(2);
    cbc[3] = new ParallelPeriodicFace<double,Dim,M,Cell>(3);
    fdi << "++++++++BConds object cbc begin++++++++" << endl;
    fdi << cbc;
    fdi << "++++++++BConds object cbc end++++++++++" << endl;

    // Cell-centered test Field's:
    Field<double,Dim,M,Cell> cA(layout,GuardCellSizes<Dim>(2),cbc);
    Field<double,Dim,M,Cell> cB(layout);

    // Set Vert-centered boundary conditions.
    BConds<double,Dim,M,Vert> vbc;
    vbc[0] = new NegReflectFace<double,Dim,M,Vert>(0);
    vbc[1] = new ZeroFace<double,Dim,M,Vert>(1);
    vbc[2] = new ParallelPeriodicFace<double,Dim,M,Vert>(2);
    vbc[3] = new ParallelPeriodicFace<double,Dim,M,Vert>(3);
    // Vert-centered test Field's:
    Field<double,Dim,M,Vert> vA(layout,GuardCellSizes<Dim>(2),vbc);
    Field<double,Dim,M,Vert> vB(layout);

    // Assign reference values:
    int i,j;
    unsigned counter=0;
    double value;
    for (j=0; j<5; j++) {
        for (i=0; i<5; i++) {
            value = counter++;
            assign(cA[i][j], value);
            assign(vA[i][j], value);
        }
    }

    // Print reference values, then assign values ofsetting across boundaries
    // and print results, Cell-centered case:
    setFormat(5,3);
    fdi << "++++++++++cA+++++++++++" << endl ;
    fp2(cA);
    cB[I][J] = cA[I-2][J-2];
    fdi << "++++++++++cB+++++++++++" << endl ;
    fp2(cB);

    // Print reference values, then assign values ofsetting across boundaries
    // and print results, Vert-centered case:
    fdi << "++++++++++vA+++++++++++" << endl ;
    fp2(vA);
    vB[I][J] = vA[I-2][J-2];
    fdi << "++++++++++vB+++++++++++" << endl ;
    fp2(vB);

    // Componentwise specification of BC's for a Field<Vektor>
    // Set Cell-centered boundary conditions.
    BConds<Vektor<double,Dim>,Dim,M,Cell> vcbc;
    vcbc[0] = new NegReflectFace<Vektor<double,Dim>,Dim,M,Cell>(0,0);
    vcbc[1] = new PosReflectFace<Vektor<double,Dim>,Dim,M,Cell>(0,1);
    vcbc[2] = new NegReflectFace<Vektor<double,Dim>,Dim,M,Cell>(1,0);
    vcbc[3] = new PosReflectFace<Vektor<double,Dim>,Dim,M,Cell>(1,1);
    vcbc[4] = new ZeroFace<Vektor<double,Dim>,Dim,M,Cell>(2);
    vcbc[5] = new ZeroFace<Vektor<double,Dim>,Dim,M,Cell>(3);
    // Cell-centered test Field's:
    Field<Vektor<double,Dim>,Dim,M,Cell> vcA(layout,GuardCellSizes<Dim>(2),vcbc);
    Field<Vektor<double,Dim>,Dim,M,Cell> vcB(layout);
    // Assign reference values:
    counter=0;
    for (j=0; j<5; j++) {
        for (i=0; i<5; i++) {
            value = counter++;
            assign(vcA[i][j], (Vektor<double,Dim>)value);
        }
    }
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, Cell-centered case:
    setFormat(2,2);
    fdi << "++++++++++vcA+++++++++++" << endl ;
    fp2(vcA);
    vcB[I][J] = vcA[I-2][J];
    fdi << "++++++++++vcB+++++++++++" << endl ;
    fp2(vcB);

    // Componentwise specification of BC's for a Field<Vektor>
    // Set CartesianCentering-centered boundary conditions.
    typedef CommonCartesianCenterings<Dim,2U>::vectorFace vFace;

    // For clarity, construct a mesh. Here, 5 is taken as the number of verts:
    M mesh(I,J);
    CenteredFieldLayout<Dim,M,vFace> layoutVFace(mesh);

    BConds<Vektor<double,Dim>,Dim,M,vFace> vfbc;
    vfbc[0] = new NegReflectFace<Vektor<double,Dim>,Dim,M,vFace>(0,0);
    vfbc[1] = new PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(0,1);
    vfbc[2] = new NegReflectFace<Vektor<double,Dim>,Dim,M,vFace>(1,0);
    vfbc[3] = new PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(1,1);
    if (Ippl::getNodes() == 1) {
        vfbc[4] = new PeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(2,0);
        vfbc[5] = new PeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(2,1);
        vfbc[6] = new PeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(3,0);
        vfbc[7] = new PeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(3,1);
    }
    else {
        vfbc[4] = new ParallelPeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(2,0);
        vfbc[5] = new ParallelPeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(2,1);
        vfbc[6] = new ParallelPeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(3,0);
        vfbc[7] = new ParallelPeriodicFace<Vektor<double,Dim>,Dim,M,vFace>(3,1);
    }
    // vFace-centered test Field's:
    Field<Vektor<double,Dim>,Dim,M,vFace>
        vfA(layoutVFace,GuardCellSizes<Dim>(2),vfbc);
    Field<Vektor<double,Dim>,Dim,M,vFace> vfB(layoutVFace);
    // Assign red-flag values for to make inaccessible vector components visible:
    vfA = 9.99;
    vfB = 9.99;
    Index Iverts(5);
    Index Jverts(5);
    Index Icells(4);
    Index Jcells(4);
    assign(vfA[Iverts][Jcells](0), Iverts + Jcells*5.0);
    assign(vfA[Icells][Jverts](1), Icells + Jverts*5.0);
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, vFace-centered case:
    setFormat(2,2);

    // Set up for 3D Field's:
    const unsigned Dim3 = 3;
    Index K(5);
    FieldLayout<Dim3> layout3(I,J,K);
    typedef UniformCartesian<Dim3> M3;
    // Componentwise specification of BC's for a Cell-centered Field<SymTenzor>
    // Set boundary conditions, positive reflecting on diagonal, negative
    // reflecting on off-diagonal:
    typedef SymTenzor<double,Dim3> ST;
    BConds<ST,Dim3,M3,Cell> tcc;
    // Face 0
    tcc[0]  = new PosReflectFace<ST,Dim3,M3,Cell>(0,0,0);
    tcc[1]  = new PosReflectFace<ST,Dim3,M3,Cell>(0,1,1);
    tcc[2]  = new PosReflectFace<ST,Dim3,M3,Cell>(0,2,2);
    tcc[3]  = new NegReflectFace<ST,Dim3,M3,Cell>(0,1,0);
    tcc[4]  = new NegReflectFace<ST,Dim3,M3,Cell>(0,2,0);
    tcc[5]  = new NegReflectFace<ST,Dim3,M3,Cell>(0,2,1);
    // Face 1
    tcc[6]  = new PosReflectFace<ST,Dim3,M3,Cell>(1,0,0);
    tcc[7]  = new PosReflectFace<ST,Dim3,M3,Cell>(1,1,1);
    tcc[8]  = new PosReflectFace<ST,Dim3,M3,Cell>(1,2,2);
    tcc[9]  = new NegReflectFace<ST,Dim3,M3,Cell>(1,1,0);
    tcc[10] = new NegReflectFace<ST,Dim3,M3,Cell>(1,2,0);
    tcc[11] = new NegReflectFace<ST,Dim3,M3,Cell>(1,2,1);
    // Face 2
    tcc[12] = new PosReflectFace<ST,Dim3,M3,Cell>(2,0,0);
    tcc[13] = new PosReflectFace<ST,Dim3,M3,Cell>(2,1,1);
    tcc[14] = new PosReflectFace<ST,Dim3,M3,Cell>(2,2,2);
    tcc[15] = new NegReflectFace<ST,Dim3,M3,Cell>(2,1,0);
    tcc[16] = new NegReflectFace<ST,Dim3,M3,Cell>(2,2,0);
    tcc[17] = new NegReflectFace<ST,Dim3,M3,Cell>(2,2,1);
    // Face 3
    tcc[18] = new PosReflectFace<ST,Dim3,M3,Cell>(3,0,0);
    tcc[19] = new PosReflectFace<ST,Dim3,M3,Cell>(3,1,1);
    tcc[20] = new PosReflectFace<ST,Dim3,M3,Cell>(3,2,2);
    tcc[21] = new NegReflectFace<ST,Dim3,M3,Cell>(3,1,0);
    tcc[22] = new NegReflectFace<ST,Dim3,M3,Cell>(3,2,0);
    tcc[24] = new NegReflectFace<ST,Dim3,M3,Cell>(3,2,1);

    // Cell-centered test Field's:
    Field<ST,Dim3,M3,Cell> sA(layout3,GuardCellSizes<Dim3>(2),tcc);
    Field<ST,Dim3,M3,Cell> sB(layout3);
    // Assign reference values:
    sA[I][J][K] = (I + J + K);
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, Cell-centered case:
    setFormat(1,2);
    fdi << "++++++++++sA+++++++++++" << endl ;
    fp3(sA);
    sB[I][J][K] = sA[I-2][J][K];
    fdi << "++++++++++sB+++++++++++" << endl ;
    fp3(sB);

    // Componentwise specification of BC's for a Field<Vektor> Set
    // CartesianCenting-centered boundary conditions.
    // TJW 12/16/97: this differs from earlier one in that it uses
    // NegFeflectZeroFace instead of NegReflectFace. This tests the new
    // NegReflectZeroFace BC, which sets last *physical* layers of vert-centered
    // quantities/components to zero. Another difference from the earlier one: I
    // put in GC and BC on the "B" and (new) "C" Field's, so that BC are applied
    // to the results of the stencil ops.
    BConds<Vektor<double,Dim>,Dim,M,vFace> vfbcz;
    vfbcz[0] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,vFace>(0,0);
    vfbcz[1] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(0,1);
    vfbcz[2] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,vFace>(1,0);
    vfbcz[3] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(1,1);
    vfbcz[4] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(2,0);
    vfbcz[5] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,vFace>(2,1);
    vfbcz[6] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,vFace>(3,0);
    vfbcz[7] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,vFace>(3,1);
    // vFace-centered test Field's:
    Field<Vektor<double,Dim>,Dim,M,vFace>
        vfzA(layout,GuardCellSizes<Dim>(2),vfbcz);
    Field<Vektor<double,Dim>,Dim,M,vFace>
        vfzB(layout,GuardCellSizes<Dim>(2),vfbcz);
    Field<Vektor<double,Dim>,Dim,M,vFace>
        vfzC(layout,GuardCellSizes<Dim>(2),vfbcz);
    // Assign reference values:
    counter=0;
    for (j=0; j<5; j++) {
        for (i=0; i<5; i++) {
            value = counter++;
            assign(vfzA[i][j], (Vektor<double,Dim>)value);
        }
    }
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, vfzace-centered case:
    setFormat(2,2);
    fdi << "++++++++++vfzA+++++++++++" << endl ;
    fp2(vfzA);
    vfzB[I][J] = vfzA[I-2][J];
    fdi << "++++++++++vfzB+++++++++++" << endl ;
    fp2(vfzB);
    vfzC[I][J] = vfzA[I-2][J-2];
    fdi << "++++++++++vfzC+++++++++++" << endl ;
    fp2(vfzC);

    // Test NegReflectZeroFace BC with Vert centering (componentwise BC):
    BConds<Vektor<double,Dim>,Dim,M,Vert> vbcz;
    vbcz[0] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,Vert>(0,0);
    vbcz[1] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,Vert>(0,1);
    vbcz[2] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,Vert>(1,0);
    vbcz[3] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,Vert>(1,1);
    vbcz[4] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,Vert>(2,0);
    vbcz[5] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,Vert>(2,1);
    vbcz[6] = new        PosReflectFace<Vektor<double,Dim>,Dim,M,Vert>(3,0);
    vbcz[7] = new NegReflectAndZeroFace<Vektor<double,Dim>,Dim,M,Vert>(3,1);
    // vFace-centered test Field's:
    Field<Vektor<double,Dim>,Dim,M,Vert> vzA(layout,GuardCellSizes<Dim>(2),vbcz);
    Field<Vektor<double,Dim>,Dim,M,Vert> vzB(layout,GuardCellSizes<Dim>(2),vbcz);
    Field<Vektor<double,Dim>,Dim,M,Vert> vzC(layout,GuardCellSizes<Dim>(2),vbcz);
    // Assign reference values:
    counter=0;
    for (j=0; j<5; j++) {
        for (i=0; i<5; i++) {
            value = counter++;
            assign(vzA[i][j], (Vektor<double,Dim>)value);
        }
    }
    // Print reference values, then assign values ofsetting across boundaries
    // and print results, vert-centered case:
    setFormat(2,2);
    fdi << "++++++++++vzA+++++++++++" << endl ;
    fp2(vzA);
    vzB[I][J] = vzA[I-2][J];
    fdi << "++++++++++vzB+++++++++++" << endl ;
    fp2(vzB);
    vzC[I][J] = vzA[I-2][J-2];
    fdi << "++++++++++vzC+++++++++++" << endl ;
    fp2(vzC);

    fdi << endl ; // Needed to flush output to file

    // Write out "by hand" into another file what the previous field-printing
    // functions should have produced; this will be compared with what they
    // actually did produce:
    hardCodedOutput("text.correct.TestBC");

    // Compare the two files by mocking up the Unix "diff" command:
    bool passed = thediff("text.test.TestBC",
                          "text.correct.TestBC");

    EXPECT_TRUE(passed);
}

namespace {
//-----------------------------------------------------------------------------
// Mock up the Unix "diff" utility to compare two files:
//-----------------------------------------------------------------------------
bool thediff(std::string filename1, std::string filename2)
{
    bool same = true;
    char ch1, ch2;
    std::ifstream file1(filename1);
    std::ifstream file2(filename2);
    while (file1.get(ch1)) {          // Read file 1 char-by-char until eof
        if (file2.get(ch2)) {           // Read equivalent char from file 2
            if (ch1 != ch2) same = false; // If they're different,files are different
        }
        else {
            same = false;                 // If file 2 ends before file 1, different
        }
    }
    return same;
}

//-----------------------------------------------------------------------------
void hardCodedOutput(std::string filename)
{
    std::ofstream of(filename);
    of << "++++++++BConds object cbc begin++++++++" << std::endl;
    of << "BConds:(" << std::endl;
    of << "NegReflectFace, Face=0 , " << std::endl;
    of << "ZeroFace, Face=1 , " << std::endl;
    of << "ParallelPeriodicFace, Face=2 , " << std::endl;
    of << "ParallelPeriodicFace, Face=3" << std::endl;
    of << ")" << std::endl;
    of << "" << std::endl;
    of << "++++++++BConds object cbc end++++++++++" << std::endl;
    of << "++++++++++cA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 9.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "1.000e+01 1.100e+01 1.200e+01 1.300e+01 1.400e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "1.500e+01 1.600e+01 1.700e+01 1.800e+01 1.900e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "2.000e+01 2.100e+01 2.200e+01 2.300e+01 2.400e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "++++++++++cB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "-1.600e+01 -1.500e+01 1.500e+01 1.600e+01 1.700e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "-2.100e+01 -2.000e+01 2.000e+01 2.100e+01 2.200e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "-1.000e+00 0.000e+00 0.000e+00 1.000e+00 2.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "-6.000e+00 -5.000e+00 5.000e+00 6.000e+00 7.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "-1.100e+01 -1.000e+01 1.000e+01 1.100e+01 1.200e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "++++++++++vA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 9.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "1.000e+01 1.100e+01 1.200e+01 1.300e+01 1.400e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "1.500e+01 1.600e+01 1.700e+01 1.800e+01 1.900e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 4.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "++++++++++vB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "-1.200e+01 -1.100e+01 1.000e+01 1.100e+01 1.200e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "-1.700e+01 -1.600e+01 1.500e+01 1.600e+01 1.700e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "-2.000e+00 -1.000e+00 0.000e+00 1.000e+00 2.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "-7.000e+00 -6.000e+00 5.000e+00 6.000e+00 7.000e+00" << std::endl;
    of << "" << std::endl << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "-1.200e+01 -1.100e+01 1.000e+01 1.100e+01 1.200e+01" << std::endl;
    of << "" << std::endl << std::endl;
    of << "++++++++++vcA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
    of << "( 2.00e+00 , 2.00e+00 ) ( 3.00e+00 , 3.00e+00 )" << std::endl;
    of << "( 4.00e+00 , 4.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 7.00e+00 , 7.00e+00 ) ( 8.00e+00 , 8.00e+00 )" << std::endl;
    of << "( 9.00e+00 , 9.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 1.20e+01 , 1.20e+01 ) ( 1.30e+01 , 1.30e+01 )" << std::endl;
    of << "( 1.40e+01 , 1.40e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 1.70e+01 , 1.70e+01 ) ( 1.80e+01 , 1.80e+01 )" << std::endl;
    of << "( 1.90e+01 , 1.90e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( 2.00e+01 , 2.00e+01 ) ( 2.10e+01 , 2.10e+01 )" << std::endl;
    of << "( 2.20e+01 , 2.20e+01 ) ( 2.30e+01 , 2.30e+01 )" << std::endl;
    of << "( 2.40e+01 , 2.40e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vcB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( -1.00e+00 , 1.00e+00 ) ( 0.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
    of << "( 2.00e+00 , 2.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( -6.00e+00 , 6.00e+00 ) ( -5.00e+00 , 5.00e+00 )" << std::endl;
    of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 7.00e+00 , 7.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( -1.10e+01 , 1.10e+01 ) ( -1.00e+01 , 1.00e+01 )" << std::endl;
    of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 1.20e+01 , 1.20e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( -1.60e+01 , 1.60e+01 ) ( -1.50e+01 , 1.50e+01 )" << std::endl;
    of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 1.70e+01 , 1.70e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( -2.10e+01 , 2.10e+01 ) ( -2.00e+01 , 2.00e+01 )" << std::endl;
    of << "( 2.00e+01 , 2.00e+01 ) ( 2.10e+01 , 2.10e+01 )" << std::endl;
    of << "( 2.20e+01 , 2.20e+01 ) " << std::endl;
    of << "" << std::endl;
    /* tjw: omit this test until componentwise ParallelPeriodicFace bug is fixed:
       of << "++++++++++vfA+++++++++++" << std::endl;
       of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
       of << "--------------------------------------------------J = 0" << std::endl;
       of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
       of << "( 2.00e+00 , 2.00e+00 ) ( 3.00e+00 , 3.00e+00 )" << std::endl;
       of << "( 4.00e+00 , 3.00e+00 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 1" << std::endl;
       of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
       of << "( 7.00e+00 , 7.00e+00 ) ( 8.00e+00 , 8.00e+00 )" << std::endl;
       of << "( 9.00e+00 , 8.00e+00 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 2" << std::endl;
       of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
       of << "( 1.20e+01 , 1.20e+01 ) ( 1.30e+01 , 1.30e+01 )" << std::endl;
       of << "( 1.40e+01 , 1.30e+01 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 3" << std::endl;
       of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
       of << "( 1.70e+01 , 1.70e+01 ) ( 1.80e+01 , 1.80e+01 )" << std::endl;
       of << "( 1.90e+01 , 1.80e+01 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 4" << std::endl;
       of << "( 9.99e+00 , 0.00e+00 ) ( 9.99e+00 , 1.00e+00 )" << std::endl;
       of << "( 9.99e+00 , 2.00e+00 ) ( 9.99e+00 , 3.00e+00 )" << std::endl;
       of << "( 9.99e+00 , 3.00e+00 ) " << std::endl;
       of << "" << std::endl;
       of << "++++++++++vfB+++++++++++" << std::endl;
       of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
       of << "--------------------------------------------------J = 0" << std::endl;
       of << "( -2.00e+00 , 1.00e+00 ) ( -1.00e+00 , 0.00e+00 )" << std::endl;
       of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
       of << "( 2.00e+00 , 2.00e+00 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 1" << std::endl;
       of << "( -7.00e+00 , 6.00e+00 ) ( -6.00e+00 , 5.00e+00 )" << std::endl;
       of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
       of << "( 7.00e+00 , 7.00e+00 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 2" << std::endl;
       of << "( -1.20e+01 , 1.10e+01 ) ( -1.10e+01 , 1.00e+01 )" << std::endl;
       of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
       of << "( 1.20e+01 , 1.20e+01 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 3" << std::endl;
       of << "( -1.70e+01 , 1.60e+01 ) ( -1.60e+01 , 1.50e+01 )" << std::endl;
       of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
       of << "( 1.70e+01 , 1.70e+01 ) " << std::endl;
       of << "" << std::endl;
       of << "--------------------------------------------------J = 4" << std::endl;
       of << "( -9.99e+00 , 1.00e+00 ) ( -9.99e+00 , 0.00e+00 )" << std::endl;
       of << "( 9.99e+00 , 0.00e+00 ) ( 9.99e+00 , 1.00e+00 )" << std::endl;
       of << "( 9.99e+00 , 2.00e+00 ) " << std::endl;
       of << "" << std::endl;
       tjw: omit this test until componentwise ParallelPeriodicFace bug is fixed. */
    of << "++++++++++sA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "==================================================K = 0" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 1" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 2" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 3" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 4" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)";
    of << "(1.10e+01 , 1.10e+01 , 1.10e+01)" << std::endl;
    of << "(1.20e+01 , 1.20e+01 , 1.20e+01)";
    of << "(1.20e+01 , 1.20e+01 , 1.20e+01)";
    of << "(1.20e+01 , 1.20e+01 , 1.20e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "++++++++++sB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "==================================================K = 0" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(1.00e+00 , -1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , 1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , -1.00e+00 , 1.00e+00)" << std::endl;
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)" << std::endl;
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)";
    of << "(0.00e+00 , 0.00e+00 , 0.00e+00)" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(2.00e+00 , -2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , 2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , -2.00e+00 , 2.00e+00)" << std::endl;
    of << "(1.00e+00 , -1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , 1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , -1.00e+00 , 1.00e+00)" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(2.00e+00 , -2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , 2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , -2.00e+00 , 2.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 1" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(2.00e+00 , -2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , 2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , -2.00e+00 , 2.00e+00)" << std::endl;
    of << "(1.00e+00 , -1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , 1.00e+00 , -1.00e+00)";
    of << "(-1.00e+00 , -1.00e+00 , 1.00e+00)" << std::endl;
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)";
    of << "(1.00e+00 , 1.00e+00 , 1.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(2.00e+00 , -2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , 2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , -2.00e+00 , 2.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 2" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(2.00e+00 , -2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , 2.00e+00 , -2.00e+00)";
    of << "(-2.00e+00 , -2.00e+00 , 2.00e+00)" << std::endl;
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)";
    of << "(2.00e+00 , 2.00e+00 , 2.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(7.00e+00 , -7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , 7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , -7.00e+00 , 7.00e+00)" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 3" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(3.00e+00 , -3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , 3.00e+00 , -3.00e+00)";
    of << "(-3.00e+00 , -3.00e+00 , 3.00e+00)" << std::endl;
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)";
    of << "(3.00e+00 , 3.00e+00 , 3.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(7.00e+00 , -7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , 7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , -7.00e+00 , 7.00e+00)" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(8.00e+00 , -8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , 8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , -8.00e+00 , 8.00e+00)" << std::endl;
    of << "(7.00e+00 , -7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , 7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , -7.00e+00 , 7.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "==================================================K = 4" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(4.00e+00 , -4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , 4.00e+00 , -4.00e+00)";
    of << "(-4.00e+00 , -4.00e+00 , 4.00e+00)" << std::endl;
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)";
    of << "(4.00e+00 , 4.00e+00 , 4.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(5.00e+00 , -5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , 5.00e+00 , -5.00e+00)";
    of << "(-5.00e+00 , -5.00e+00 , 5.00e+00)" << std::endl;
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)";
    of << "(5.00e+00 , 5.00e+00 , 5.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "(7.00e+00 , -7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , 7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , -7.00e+00 , 7.00e+00)" << std::endl;
    of << "(6.00e+00 , -6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , 6.00e+00 , -6.00e+00)";
    of << "(-6.00e+00 , -6.00e+00 , 6.00e+00)" << std::endl;
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)";
    of << "(6.00e+00 , 6.00e+00 , 6.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "(8.00e+00 , -8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , 8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , -8.00e+00 , 8.00e+00)" << std::endl;
    of << "(7.00e+00 , -7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , 7.00e+00 , -7.00e+00)";
    of << "(-7.00e+00 , -7.00e+00 , 7.00e+00)" << std::endl;
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)";
    of << "(7.00e+00 , 7.00e+00 , 7.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "(9.00e+00 , -9.00e+00 , -9.00e+00)";
    of << "(-9.00e+00 , 9.00e+00 , -9.00e+00)";
    of << "(-9.00e+00 , -9.00e+00 , 9.00e+00)" << std::endl;
    of << "(8.00e+00 , -8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , 8.00e+00 , -8.00e+00)";
    of << "(-8.00e+00 , -8.00e+00 , 8.00e+00)" << std::endl;
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)";
    of << "(8.00e+00 , 8.00e+00 , 8.00e+00)" << std::endl;
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)";
    of << "(9.00e+00 , 9.00e+00 , 9.00e+00)" << std::endl;
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)";
    of << "(1.00e+01 , 1.00e+01 , 1.00e+01)" << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;
    of << "++++++++++vfzA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
    of << "( 2.00e+00 , 2.00e+00 ) ( 3.00e+00 , 3.00e+00 )" << std::endl;
    of << "( -3.00e+00 , 3.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 7.00e+00 , 7.00e+00 ) ( 8.00e+00 , 8.00e+00 )" << std::endl;
    of << "( -8.00e+00 , 8.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 1.20e+01 , 1.20e+01 ) ( 1.30e+01 , 1.30e+01 )" << std::endl;
    of << "( -1.30e+01 , 1.30e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 1.70e+01 , 1.70e+01 ) ( 1.80e+01 , 1.80e+01 )" << std::endl;
    of << "( -1.80e+01 , 1.80e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( 1.50e+01 , -1.50e+01 ) ( 1.60e+01 , -1.60e+01 )" << std::endl;
    of << "( 1.70e+01 , -1.70e+01 ) ( 1.80e+01 , -1.80e+01 )" << std::endl;
    of << "( -1.80e+01 , -1.80e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vfzB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( -1.00e+00 , 1.00e+00 ) ( 0.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
    of << "( -1.00e+00 , 1.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( -6.00e+00 , 6.00e+00 ) ( -5.00e+00 , 5.00e+00 )" << std::endl;
    of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( -6.00e+00 , 6.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( -1.10e+01 , 1.10e+01 ) ( -1.00e+01 , 1.00e+01 )" << std::endl;
    of << "( 1.00e+01 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( -1.10e+01 , 1.10e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( -1.60e+01 , 1.60e+01 ) ( -1.50e+01 , 1.50e+01 )" << std::endl;
    of << "( 1.50e+01 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( -1.60e+01 , 1.60e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( -1.60e+01 , -1.60e+01 ) ( -1.50e+01 , -1.50e+01 )" << std::endl;
    of << "( 1.50e+01 , -1.50e+01 ) ( 1.60e+01 , -1.60e+01 )" << std::endl;
    of << "( -1.60e+01 , -1.60e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vfzC+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( -6.00e+00 , -6.00e+00 ) ( -5.00e+00 , -5.00e+00 )" << std::endl;
    of << "( 5.00e+00 , -5.00e+00 ) ( 6.00e+00 , -6.00e+00 )" << std::endl;
    of << "( -6.00e+00 , -6.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( -1.00e+00 , -1.00e+00 ) ( 0.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , -1.00e+00 )" << std::endl;
    of << "( -1.00e+00 , -1.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( -1.00e+00 , 1.00e+00 ) ( 0.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 1.00e+00 )" << std::endl;
    of << "( -1.00e+00 , 1.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( -6.00e+00 , 6.00e+00 ) ( -5.00e+00 , 5.00e+00 )" << std::endl;
    of << "( 5.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( -6.00e+00 , 6.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( -6.00e+00 , -6.00e+00 ) ( -5.00e+00 , -5.00e+00 )" << std::endl;
    of << "( 5.00e+00 , -5.00e+00 ) ( 6.00e+00 , -6.00e+00 )" << std::endl;
    of << "( -6.00e+00 , -6.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vzA+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 2.00e+00 , 0.00e+00 ) ( 3.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( 0.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 7.00e+00 , 7.00e+00 ) ( 8.00e+00 , 8.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 9.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( 0.00e+00 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 1.20e+01 , 1.20e+01 ) ( 1.30e+01 , 1.30e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.40e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( 0.00e+00 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 1.70e+01 , 1.70e+01 ) ( 1.80e+01 , 1.80e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.90e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 2.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 2.20e+01 , 0.00e+00 ) ( 2.30e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vzB+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( -1.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( 0.00e+00 , 7.00e+00 ) ( -6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 7.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( 0.00e+00 , 1.20e+01 ) ( -1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.00e+01 ) ( 1.10e+01 , 1.10e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.20e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( 0.00e+00 , 1.70e+01 ) ( -1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.50e+01 ) ( 1.60e+01 , 1.60e+01 )" << std::endl;
    of << "( 0.00e+00 , 1.70e+01 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( -2.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 2.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "++++++++++vzC+++++++++++" << std::endl;
    of << "~~~~~~~~ field slice (0:4:1, 0:4:1) ~~~~~~~~" << std::endl;
    of << "--------------------------------------------------J = 0" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( -1.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 1" << std::endl;
    of << "( 0.00e+00 , -7.00e+00 ) ( -6.00e+00 , -6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , -5.00e+00 ) ( 6.00e+00 , -6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , -7.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 2" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( -1.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.00e+00 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 3" << std::endl;
    of << "( 0.00e+00 , 7.00e+00 ) ( -6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 5.00e+00 ) ( 6.00e+00 , 6.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 7.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "--------------------------------------------------J = 4" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( -1.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) ( 1.10e+01 , 0.00e+00 )" << std::endl;
    of << "( 0.00e+00 , 0.00e+00 ) " << std::endl;
    of << "" << std::endl;
    of << "" << std::endl;

    of.close();
}
}