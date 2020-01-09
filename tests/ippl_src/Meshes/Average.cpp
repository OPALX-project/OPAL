// Adapted from TestAverageVC_CV.cpp , Tim Williams 2/6/1997
// This tests the Average functions in [Uniform]Cartesian class, which do
// weighted averages between Cell and Vertex centered Field's.

#include "gtest/gtest.h"

#include "Meshes/UniformCartesian.h"

// set dimensionality and problem size
const unsigned Dim3 = 3;
const unsigned nx = 4, ny = 4, nz = 4;
constexpr double roundOffError = 1e-10;

using namespace IPPL;

TEST(Meshes, Average)
{
    GuardCellSizes<Dim3> gc(1);
    Index I(nx);
    Index J(ny);
    Index K(nz);
    FieldLayout<Dim3> layoutVert(I,J,K);
    Index Ic(nx-1);
    Index Jc(ny-1);
    Index Kc(nz-1);
    FieldLayout<Dim3> layoutCell(Ic,Jc,Kc);

    // Scalar Field, Cell-Centered
    Field<double,Dim3,UniformCartesian<Dim3>,Cell> C(layoutCell,gc);
    C = 1.0;

    // Scalar weight Field, Cell-Centered
    Field<double,Dim3,UniformCartesian<Dim3>,Cell> wC(layoutCell,gc);
    wC = 2.0;

    // Scalar Field, Vert-Centered
    Field<double,Dim3,UniformCartesian<Dim3>,Vert> V(layoutVert,gc);
    V = 1.0;

    // Scalar weight Field, Vert-Centered
    Field<double,Dim3,UniformCartesian<Dim3>,Vert> wV(layoutVert,gc);
    wV = 2.0;

    // Field's to hold weighted averages:
    Field<double,Dim3,UniformCartesian<Dim3>,Cell> avgToC(layoutCell,gc);
    Field<double,Dim3,UniformCartesian<Dim3>,Vert> avgToV(layoutVert,gc);

    // Weighted average from Cell to Vert:
    assign(avgToV, Average(C, wC, avgToV));

    // Weighted average from Vert to Cell:
    assign(avgToC, Average(V, wV, avgToC));

    // Weight from Field<Vektor,Vert> to Field<Vektor,Cell>, using scalar
    // weight field wV (Field<double,Vert>)
    // Vector Field, Cell-centered
    Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Vert>
        vV(layoutVert,gc);
    Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell>
        avgToCSW(layoutCell,gc);
    Vektor<double,Dim3> ones;
    for (unsigned int d=0; d < Dim3; d++) ones(d) = 1.0;
    vV = ones;
    assign(avgToCSW, Average(vV, wV, avgToCSW));

    // Check results:
    EXPECT_NEAR(sum(avgToV), nx*ny*nz, roundOffError);
    EXPECT_NEAR(sum(avgToC), (nx-1)*(ny-1)*(nz-1), roundOffError);

    // Following triggers two IPPL bugs:
    // 1) unsigned*Vektor<double,Dim> doesn't work, have to cast the
    //    unsigned to a double
    // 2) sum(Field<Vektor,...>,...>) doesn't work.
    //  if (sum(avgToCSW) != (nx-1)*(ny-1)*(nz-1)*ones) passed = false;
    // WORKAROUND:
    Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell>::iterator fi;
    Vektor<double,Dim3> thesum = 0.0;
    Vektor<double,Dim3> globalSum = 0.0;
    for (fi = avgToCSW.begin(); fi != avgToCSW.end(); ++fi) thesum += *fi;
    // This gets the global sum to every PE's copy when multiprocessing:
    globalSum = thesum;
    reduce(globalSum, globalSum, OpAddAssign());
    EXPECT_TRUE(globalSum == ((double)(nx-1)*(ny-1)*(nz-1))*ones);

    // Test one of the unweighted (2-argument) Average() functions:
    // Average from Field<Vektor,Vert> to Field<Vektor,Cell>
    // Vector Field, Cell-centered
    Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Vert>
        vVu(layoutVert,gc);
    Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Cell>
        avgToCSu(layoutCell,gc);
    vVu = ones;
    assign(avgToCSu, Average(vV, avgToCSu));
    // Check results:
    EXPECT_TRUE(sum(avgToCSu) == (nx-1)*(ny-1)*(nz-1));
}