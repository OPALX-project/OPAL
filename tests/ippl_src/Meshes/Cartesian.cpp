#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Meshes/Cartesian.h"

#include "FieldLayout/CenteredFieldLayout.h"

constexpr unsigned int DIM = 3;
constexpr double roundOffError = 1e-10;

TEST(Meshes, Cartesian)
{
    OpalTestUtilities::SilenceTest silencer;

    const unsigned D = DIM; // Hardwire dimensionality
    const unsigned nv = 6;  // Hardwire number of vertices in every direction
    unsigned vnodes = 4;    // Hardwire 4 vnodes

    // Sizes:
    unsigned nverts[D], ncells[D];
    unsigned totverts=1, totcells=1;
    unsigned int d;

    for (d=0; d<D; d++) {
        ncells[d] = nv - 1;
        nverts[d] = nv;
        totcells *= ncells[d];
        totverts *= nverts[d];
    }
    NDIndex<D> verts, cells;
    for (d=0; d<D; d++) {
        verts[d] = Index(nverts[d]);
        cells[d] = Index(ncells[d]);
    }

    //---------------------------------------------------------------------------
    // Construct some CenteredFieldLayout's and Field's to be used below:

    // Create cartesian mesh object:
    typedef Cartesian<D,double> M;

    double* delX[D];

    for (d=0; d<D; d++)
        delX[d] = new double[nverts[d]];

    Vektor<double,D> origin;
    for (d=0; d<D; d++)
        origin(d) = d + 1.0;

    // Assign nonuniform mesh-spacing values to each component (linear ramps):
    for (d=0; d<D; d++) {
        double multipplier = (d + 1)*1.0;
        for (unsigned int vert=0; vert < nverts[d]; vert++) {
            (delX[d])[vert] = multipplier*(1 + vert);
        }
    }

    // Mesh boundary conditions:
    MeshBC_E mbc[2*D];
    for (unsigned b=0; b < (2*D); b++)
        mbc[b] = Reflective;

    // Test constructing mesh, and then setting spacing, origin, BC's
    M mesh(verts);
    mesh.set_origin(origin);
    mesh.set_meshSpacing(delX);
    mesh.set_MeshBC(mbc);

    // Clean up mesh spacing arrays
    for (d=0; d<D; d++)
        delete [] delX[d];

    // ada have to cross check Div() fails without this
    mesh.storeSpacingFields();

    // Construct CenteredFieldLayout's using this for Vert and Cell centering:
    e_dim_tag edt[D];
    for (d=0; d<D; d++)
        edt[d] = PARALLEL;
    CenteredFieldLayout<D,M,Cell> cl(mesh, edt, vnodes);
    CenteredFieldLayout<D,M,Vert> vl(mesh, edt, vnodes);

    // Use 1 guard layer in all Field's:
    GuardCellSizes<D> gc(1);

    // Vectors:
    BConds<Vektor<double,D>,D,M,Vert> vvbc;
    BConds<Vektor<double,D>,D,M,Cell> vcbc;

    // Scalars:
    BConds<double,D,M,Cell> scbc;

    // Symmetric tensors:
    BConds<SymTenzor<double,D>,D,M,Cell> stcbc;

    // Tensors:
    BConds<Tenzor<double,D>,D,M,Cell> tcbc;

    // Use linear negative reflecting conditions:
    for (unsigned int face=0; face<2*D; face++) {
        vvbc[face]  = new NegReflectFace<Vektor<double,D>,D,M,Vert>(face);
        vcbc[face]  = new NegReflectFace<Vektor<double,D>,D,M,Cell>(face);
        scbc[face]  = new NegReflectFace<double,D,M,Cell>(face);
        stcbc[face] = new NegReflectFace<SymTenzor<double,D>,D,M,Cell>(face);
        tcbc[face] =  new NegReflectFace<Tenzor<double,D>,D,M,Cell>(face);
    }

    // Now use all this to construct some Field's:
    Field<Vektor<double,D>,D,M,Vert> vectorVert(mesh, vl, gc, vvbc);
    Field<Vektor<double,D>,D,M,Cell> vectorCell(mesh, cl, gc, vcbc);
    Field<SymTenzor<double,D>,D,M,Cell> symtCell(mesh, cl, gc, stcbc);
    Field<Tenzor<double,D>,D,M,Cell> tensorCell(mesh, cl, gc, tcbc);
    Field<double,D,M,Cell> scalarCell(mesh, cl, gc, scbc);

    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Try out Divergence Vektor/Vert -> Scalar/Cell:
    // Assign values into the vert-centered Field<Vektor>:
    assign(vectorVert, mesh.getVertexPositionField(vectorVert));
    scalarCell = Div(vectorVert, scalarCell);
    // The value should be 3.0 for all elements; test this:
    EXPECT_NEAR(abs(sum(scalarCell)/totcells), 1.0*D, roundOffError);
    //---------------------------------------------------------------------------

    // --------------------------------------------------------------------------
    // Try out Gradient Scalar/Cell -> Vektor/Vert:

    // Use mesh object and vectorVert and scalarCell Field's constructed above.
    vectorCell = mesh.getCellPositionField(vectorCell);
    vectorCell -= mesh.get_origin();
    // Assign positive-sloping linear ramp values into the cell-centered
    // Field<scalar>:
    scalarCell = 0.0;
    for (d=0; d<D; d++) scalarCell[cells] += vectorCell[cells](d);
    // Now take the gradient:
    vectorVert = Grad(scalarCell, vectorVert);
    // The value should be (1.0,1.0,1.0) for all elements one at least one
    // removed from the last-physical-layer elements. Last-physical-layer
    // elements will be different because the BC available in IPPL don't really
    // do the kind of linear extrapolation appropriate for the needs here:
    Vektor<double,D> unit; for (d=0; d<D; d++) unit[d] = 1.0;
    Vektor<double,D> sumVectorVert;
    // Use temporary, smaller BareField as a reduced-by-two vector Field to hold
    // only the boundary-exclusive elements (needed because of limitations of
    // IPPL reductions ops):
    NDIndex<D> bev;
    for (d=0; d<D; d++) bev[d] = Index(1,nverts[d]-2,1);
    FieldLayout<D> templayout(bev);
    BareField<Vektor<double,D>,D> temp(templayout);
    temp[bev] = vectorVert[bev];
    sumVectorVert = sum(temp);
    unsigned totred=1; for (d=0; d<D; d++) totred *= nverts[d] - 2;
    sumVectorVert /= totred;
    Vektor<double,D> diffVectorVert;
    diffVectorVert = sumVectorVert - unit;
    double magDiffVectorVert = 0.0;
    for (d=0; d<D; d++) magDiffVectorVert += diffVectorVert(d)*diffVectorVert(d);
    magDiffVectorVert = sqrt(magDiffVectorVert);
    EXPECT_NEAR(abs(magDiffVectorVert), 0, roundOffError);
    //---------------------------------------------------------------------------

    // --------------------------------------------------------------------------
    // Try out Gradient Scalar/Cell -> Vektor/Cell:

    // Use mesh object and vectorVert and scalarCell Field's constructed above.
    vectorCell = mesh.getCellPositionField(vectorCell);
    vectorCell -= mesh.get_origin();
    // Assign positive-sloping linear ramp values into the cell-centered
    // Field<scalar>:
    scalarCell = 0.0;
    for (d=0; d<D; d++) scalarCell[cells] += vectorCell[cells](d);
    // Now take the gradient:
    vectorCell = Grad(scalarCell, vectorCell);
    // The value should be (1.0,1.0,1.0) for all elements one at least one
    // removed from the last-physical-layer elements. Last-physical-layer
    // elements will be different because the BC available in IPPL don't really
    // do the kind of linear extrapolation appropriate for the needs here:
    for (d=0; d<D; d++) unit[d] = 1.0;
    Vektor<double,D> sumVectorCell;
    // Use temporary, smaller BareField as a reduced-by-two vector Field to hold
    // only the boundary-exclusive elements (needed because of limitations of
    // IPPL reductions ops):
    NDIndex<D> bec;
    for (d=0; d<D; d++) bec[d] = Index(1,ncells[d]-2,1);
    FieldLayout<D> templayout2(bec);
    BareField<Vektor<double,D>,D> temp2(templayout);
    temp2[bec] = vectorCell[bec];
    sumVectorCell = sum(temp2);
    unsigned totredc=1; for (d=0; d<D; d++) totredc *= ncells[d] - 2;
    sumVectorCell /= totredc;
    Vektor<double,D> diffVectorCell;
    diffVectorCell = sumVectorCell - unit;
    double magDiffVectorCell = 0.0;
    for (d=0; d<D; d++) magDiffVectorCell += diffVectorCell(d)*diffVectorCell(d);
    magDiffVectorCell = sqrt(magDiffVectorCell);
    EXPECT_NEAR(abs(magDiffVectorCell), 0, roundOffError);
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // Try out Divergence SymTenzor/Cell -> Vektor/Vert:

    // Use CenteredFieldLayout's from above object to construct SymTenzor Field:
    // Assign values into the cell-centered Field<SymTenzor>; use values from
    // cell-centered scalar Field scalarCell set up above:
    SymTenzor<double,D> unitSymTenzor = 1.0;
    symtCell = unitSymTenzor*scalarCell;
    // Now take the divergence:
    vectorVert = Div(symtCell, vectorVert);
    // The value should be (D,D,D,....) for all elements; test this:
    // Use temporary, smaller BareField as a reduced-by-two symtensor Field to
    // hold only the boundary-exclusive elements (needed because of limitations
    // of IPPL reductions ops):
    temp[bev] = vectorVert[bev];
    sumVectorVert = sum(temp);
    sumVectorVert /= totred;
    Vektor<double,D> deesVector; for (d=0; d<D; d++) deesVector(d) = 1.0*D;
    diffVectorVert = sumVectorVert - deesVector;
    magDiffVectorVert = 0.0;
    for (d=0; d<D; d++) magDiffVectorVert += diffVectorVert(d)*diffVectorVert(d);
    magDiffVectorVert = sqrt(magDiffVectorVert);
    EXPECT_NEAR(abs(magDiffVectorCell), 0, roundOffError);
    //---------------------------------------------------------------------------

    // --------------------------------------------------------------------------
    // Try out Gradient Vektor/Vert -> Tenzor/Cell:

    // Set up input values in Vektor/Vert field:
    vectorVert = mesh.getVertexPositionField(vectorVert);
    // Now take the gradient:
    tensorCell = Grad(vectorVert, tensorCell);
    // Since this is the gradient of the position vector (x*x_hat + y* y_hat +
    // z*z_hat), the result should be the identity tensor (NRL Plasma Formulary
    // Vector Identities section):
    Tenzor<double,D> identityTensor = 0.0;
    for (d=0; d<D; d++) identityTensor(d,d) = 1.0;
    Tenzor<double,D> sumTensorCell = sum(tensorCell);
    sumTensorCell /= totcells;
    Tenzor<double,D> diffTensorCell;
    diffTensorCell = sumTensorCell - identityTensor;
    double magDiffTensorCell = 0.0;
    for (d=0; d<D; d++) {
        for (unsigned int d2=0; d2<D; d2++) {
            magDiffTensorCell += diffTensorCell(d,d2)*diffTensorCell(d,d2);
        }
    }
    EXPECT_NEAR(abs(magDiffVectorCell), 0, roundOffError);
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------

    /* THIS TEST DOES NOT COMPILE

    // Test Average() functions:
    // Scalar Field, Cell-Centered
    Field<double,D,Cartesian<D>,Cell> C(mesh, cl, gc);
    C = 1.0;
    // Scalar weight Field, Cell-Centered
    Field<double,D,Cartesian<D>,Cell> wC(mesh, cl, gc);
    wC = 2.0;
    // Scalar Field, Vert-Centered
    Field<double,D,Cartesian<D>,Vert> V(mesh, vl, gc);
    V = 1.0;
    // Scalar weight Field, Vert-Centered
    Field<double,D,Cartesian<D>,Vert> wV(mesh, vl, gc);
    wV = 2.0;
    // Field's to hold weighted averages:
    Field<double,D,Cartesian<D>,Cell> avgToC(mesh, cl, gc);
    Field<double,D,Cartesian<D>,Vert> avgToV(mesh, vl, gc);

    assign(avgToV, Average(C, wC, avgToV));
    assign(avgToC, Average(V, wV, avgToC));

    // Weighted average from Cell to Vert:
    // ada  does not work assign(avgToV, Average(C, wC, avgToV));
    // Weighted average from Vert to Cell:
    // ada dones not work assign(avgToC, Average(V, wV, avgToC));
    // Check results:
    if (sum(avgToV) != totverts) {
    testmsg << "avgToV values wrong" << endl;
    testmsg << "sum(avgToV) = " << sum(avgToV) << " ; totverts = " << totverts
    << endl;
    // passed = false;
    }
    if (sum(avgToC) != totcells) {
    testmsg << "avgToC values wrong" << endl;
    testmsg << "sum(avgToC) = " << sum(avgToC) << " ; totcells = " << totcells
    << endl;
    // passed = false;
    }

    */

    //---------------------------------------------------------------------------

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //---------------------------------------------------------------------------
    // Some accessor function tests:
    double theVolume, theVolume2, theVolume3;
    NDIndex<D> ndi;
    ndi[0] = Index(2,2,1);
    ndi[1] = Index(2,2,1);
    ndi[2] = Index(2,2,1);
    theVolume = mesh.getCellVolume(ndi);
    ndi[0] = Index(0,2,1);
    ndi[1] = Index(0,2,1);
    ndi[2] = Index(0,2,1);
    theVolume2 = mesh.getCellRangeVolume(ndi);
    EXPECT_NEAR(theVolume2, (6.0*12.0*18.0), roundOffError);

    ndi[0] = Index(0,3,1);
    ndi[1] = Index(0,3,1);
    ndi[2] = Index(0,3,1);
    theVolume3 = mesh.getVertRangeVolume(ndi);
    EXPECT_NEAR(theVolume3, theVolume2, roundOffError);
    //---------------------------------------------------------------------------
    Field<double,D,Cartesian<D>,Cell> theVolumes(mesh, cl);
    mesh.getCellVolumeField(theVolumes);
    EXPECT_NEAR((sum(theVolumes)/totcells), theVolume, roundOffError);
    //---------------------------------------------------------------------------
    Vektor<double,D> v;
    v(0) = 1.5; v(1) = 4.5; v(2) = 9.5;
    ndi = mesh.getNearestVertex(v);
    // nearest vertex should be (1,1,1)
    for (unsigned int i=0; i<D; i++) {
        EXPECT_EQ((int)ndi[0].first(),  1);
        EXPECT_EQ((int)ndi[0].length(), 1);
    }
    //---------------------------------------------------------------------------
    Vektor<double,D> v1;
    v1 = mesh.getVertexPosition(ndi);
    v(0) = 2.0; v(1) = 4.0; v(2) = 12.0; // Correct value
    for (unsigned int i=0; i<D; i++) {
        EXPECT_NEAR(v1(i), v(i), roundOffError);
    }
    //---------------------------------------------------------------------------
    CenteredFieldLayout<D,Cartesian<D>,Vert>
        clVert(mesh);
    Field<Vektor<double,D>,D,Cartesian<D>,Vert>
        thePositions(clVert);
    mesh.getVertexPositionField(thePositions);
    //---------------------------------------------------------------------------
    v = mesh.getDeltaVertex(ndi);
    Vektor<double,D> vcorrect;
    vcorrect(0) = 2.0; vcorrect(1) = 4.0; vcorrect(2) = 9.0;
    for (unsigned int i=0; i<D; i++) {
        EXPECT_NEAR(vcorrect(i), v(i), roundOffError);
    }
}