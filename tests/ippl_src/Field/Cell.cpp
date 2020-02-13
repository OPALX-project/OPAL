#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Field/BareField.h"
#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Index/Index.h"
#include "Meshes/UniformCartesian.h"
#include "Utility/FieldDebug.h"

#include <iostream>

constexpr double   roundOffError = 1e-10;

namespace {
    //----------------------------------------------------------------------
    /*
      void
      initTenzorCell(Field<Tenzor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1) {
      double X[1] = { 1.0 } ;
      C1.get_mesh().set_meshSpacing(X);
      FieldLayout<1U>& layoutCell = C1.getLayout();
      const GuardCellSizes<1U>& guardCells = C1.getGuardCellSizes();
      Field<double,1U,UniformCartesian<1U>,Cell>
      Vscalar1(layoutCell, guardCells);
      Field <Tenzor<double,1U>,1U,UniformCartesian<1U>,Cell >::iterator pv;
      Field <double,1U,UniformCartesian<1U>,Cell >::iterator pf1;
      Index I = C1.getIndex(0);
      Vscalar1[I]= I;
      C1.Uncompress();
      for (pv=C1.begin(), pf1=Vscalar1.begin() ;
      pv!=C1.end() ; ++pv, ++pf1) {
      (*pv)(0,0) = *pf1;
      }
      }
    */
    //----------------------------------------------------------------------
    void
    initTenzorCell(Field<Tenzor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1) {
        double X[2] = { 1.0, 1.0 } ;
        C1.get_mesh().set_meshSpacing(X);
        FieldLayout<2U>& layoutCell = C1.getLayout();
        const GuardCellSizes<2U>& guardCells = C1.getGuardCellSizes();
        Field<double,2U,UniformCartesian<2U>,Cell>
            Vscalar1(layoutCell, guardCells),
            Vscalar2(layoutCell, guardCells);
        Field <Tenzor<double,2U>,2U,UniformCartesian<2U>,Cell >::iterator pv;
        Field <double,2U,UniformCartesian<2U>,Cell >::iterator pf1;
        Index I = C1.getIndex(0);
        Index J = C1.getIndex(1);
        Vscalar1[I][J]= I;
        C1.Uncompress();
        for (pv=C1.begin(), pf1=Vscalar1.begin() ; pv!=C1.end() ; ++pv, ++pf1) {
            (*pv)(0,0) = *pf1;
            (*pv)(0,1) = 0.0;
            (*pv)(1,0) = 0.0;
            (*pv)(1,1) = *pf1;
        }
    }
    //----------------------------------------------------------------------
    /*
      void
      initTenzorCell(Field<Tenzor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1) {
      double X[3] = { 1.0, 1.0, 1.0 } ;
      C1.get_mesh().set_meshSpacing(X);
      FieldLayout<3U>& layoutCell = C1.getLayout();
      const GuardCellSizes<3U>& guardCells = C1.getGuardCellSizes();
      Field<double,3U,UniformCartesian<3U>,Cell>
      Vscalar1(layoutCell, guardCells),
      Vscalar2(layoutCell, guardCells),
      Vscalar3(layoutCell, guardCells);
      Field <Tenzor<double,3U>,3U,UniformCartesian<3U>,Cell >::iterator pv;
      Field <double,3U,UniformCartesian<3U>,Cell >::iterator pf1, pf2, pf3;
      Index I = C1.getIndex(0);
      Index J = C1.getIndex(1);
      Index K = C1.getIndex(2);
      Vscalar1[I][J][K]= I;
      C1.Uncompress();
      for (pv=C1.begin(), pf1=Vscalar1.begin() ; pv!=C1.end() ; ++pv, ++pf1) {
      (*pv)(0,0) = *pf1;
      (*pv)(0,1) = 0.0;
      (*pv)(0,2) = 0.0;
      (*pv)(1,0) = 0.0;
      (*pv)(1,1) = *pf1;
      (*pv)(1,2) = 0.0;
      (*pv)(2,0) = 0.0;
      (*pv)(2,1) = 0.0;
      (*pv)(2,2) = *pf1;
      }
      }
    */
    //----------------------------------------------------------------------
    /*
      void
      initVectorCell(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1) {
      double X[1] = { 1.0 } ;
      C1.get_mesh().set_meshSpacing(X);
      FieldLayout<1U>& layoutCell = C1.getLayout();
      const GuardCellSizes<1U>& guardCells = C1.getGuardCellSizes();
      Field<double,1U,UniformCartesian<1U>,Cell>
      Vscalar1(layoutCell, guardCells);
      Field <Vektor<double,1U>,1U,UniformCartesian<1U>,Cell >::iterator pv;
      Field <double,1U,UniformCartesian<1U>,Cell >::iterator pf1;
      Index I = C1.getIndex(0);
      Vscalar1[I]= I;
      C1.Uncompress();
      for (pv=C1.begin(), pf1=Vscalar1.begin() ;
      pv!=C1.end() ; ++pv, ++pf1) {
      (*pv)[0] = *pf1;
      }
      }
    */
    //----------------------------------------------------------------------
    void
    initVectorCell(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1) {
        double X[2] = { 1.0, 1.0 } ;
        C1.get_mesh().set_meshSpacing(X);
        FieldLayout<2U>& layoutCell = C1.getLayout();
        const GuardCellSizes<2U>& guardCells = C1.getGuardCellSizes();
        Field<double,2U,UniformCartesian<2U>,Cell>
            Vscalar1(layoutCell, guardCells),
            Vscalar2(layoutCell, guardCells);
        Field <Vektor<double,2U>,2U,UniformCartesian<2U>,Cell >::iterator pv;
        Field <double,2U,UniformCartesian<2U>,Cell >::iterator pf1, pf2;
        Index I = C1.getIndex(0);
        Index J = C1.getIndex(1);
        Vscalar1[I][J]= I;
        Vscalar2[I][J]= J;
        C1.Uncompress();
        for (pv=C1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ;
             pv!=C1.end() ; ++pv, ++pf1, ++pf2) {
            (*pv)[0] = *pf1;
            (*pv)[1] = *pf2;
        }
    }
    //----------------------------------------------------------------------
    void
    initVectorCell(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1) {
        double X[3] = { 1.0, 1.0, 1.0 } ;
        C1.get_mesh().set_meshSpacing(X);
        FieldLayout<3U>& layoutCell = C1.getLayout();
        const GuardCellSizes<3U>& guardCells = C1.getGuardCellSizes();
        Field<double,3U,UniformCartesian<3U>,Cell>
            Vscalar1(layoutCell, guardCells),
            Vscalar2(layoutCell, guardCells),
            Vscalar3(layoutCell, guardCells);
        Field <Vektor<double,3U>,3U,UniformCartesian<3U>,Cell >::iterator pv;
        Field <double,3U,UniformCartesian<3U>,Cell >::iterator pf1, pf2, pf3;
        Index I = C1.getIndex(0);
        Index J = C1.getIndex(1);
        Index K = C1.getIndex(2);
        Vscalar1[I][J][K]= I;
        Vscalar2[I][J][K]= J;
        Vscalar3[I][J][K]= K;
        C1.Uncompress();
        for (pv=C1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ,
                 pf3=Vscalar3.begin() ;
             pv!=C1.end() ; ++pv, ++pf1, ++pf2, ++pf3) {
            (*pv)[0] = *pf1;
            (*pv)[1] = *pf2;
            (*pv)[2] = *pf3;
        }
    }
    //----------------------------------------------------------------------
    /*
      void initScalarCell(Field<double,1U,UniformCartesian<1U>,Cell>& C1) {
      Index I = C1.getIndex(0);
      C1[I]= I ;
      }
    */
    //----------------------------------------------------------------------
    void initScalarCell(Field<double,2U,UniformCartesian<2U>,Cell>& C1) {
        Index I = C1.getIndex(0);
        Index J = C1.getIndex(1);
        C1[I][J]= I ;
    }
    //----------------------------------------------------------------------
    void initScalarCell(Field<double,3U,UniformCartesian<3U>,Cell>& C1) {
        Index I = C1.getIndex(0);
        Index J = C1.getIndex(1);
        Index K = C1.getIndex(2);
        C1[I][J][K]= I ;
    }

    //----------------------------------------------------------------------
    /*
    void initVectorVert(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Vert>& V1) {
        double X[1] = { 1.0 } ;
        V1.get_mesh().set_meshSpacing(X);
        FieldLayout<1U>& layoutVert = V1.getLayout();
        const GuardCellSizes<1U>& guardCells = V1.getGuardCellSizes();
        Field<double,1U,UniformCartesian<1U>,Vert>
            Vscalar1(layoutVert, guardCells);
        Field <Vektor<double,1U>,1U,UniformCartesian<1U>,Vert >::iterator pv;
        Field <double,1U,UniformCartesian<1U>,Vert >::iterator pf1;
        Index I = V1.getIndex(0);
        Vscalar1[I]= I;
        V1.Uncompress();
        for (pv=V1.begin(), pf1=Vscalar1.begin();
             pv!=V1.end() ; ++pv, ++pf1) {
            (*pv)[0] = *pf1;
        }
    }
    //----------------------------------------------------------------------
    void initVectorVert(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Vert>& V1) {
        double X[2] = { 1.0, 1.0 } ;
        V1.get_mesh().set_meshSpacing(X);
        FieldLayout<2U>& layoutVert = V1.getLayout();
        const GuardCellSizes<2U>& guardCells = V1.getGuardCellSizes();
        Field<double,2U,UniformCartesian<2U>,Vert>
            Vscalar1(layoutVert, guardCells),
            Vscalar2(layoutVert, guardCells);
        Field <Vektor<double,2U>,2U,UniformCartesian<2U>,Vert >::iterator pv;
        Field <double,2U,UniformCartesian<2U>,Vert >::iterator pf1, pf2;
        Index I = V1.getIndex(0);
        Index J = V1.getIndex(1);
        Vscalar1[I][J]= I;
        Vscalar2[I][J]= J;
        V1.Uncompress();
        for (pv=V1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ;
             pv!=V1.end() ; ++pv, ++pf1, ++pf2) {
            (*pv)[0] = *pf1;
            (*pv)[1] = *pf2;
        }
    }
    */
    //----------------------------------------------------------------------
    void initVectorVert(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Vert>& V1) {
        double X[3] = { 1.0, 1.0, 1.0 } ;
        V1.get_mesh().set_meshSpacing(X);
        FieldLayout<3U>& layoutVert = V1.getLayout();
        const GuardCellSizes<3U>& guardCells = V1.getGuardCellSizes();
        Field<double,3U,UniformCartesian<3U>,Vert>
            Vscalar1(layoutVert, guardCells),
            Vscalar2(layoutVert, guardCells),
            Vscalar3(layoutVert, guardCells);
        Field <Vektor<double,3U>,3U,UniformCartesian<3U>,Vert >::iterator pv;
        Field <double,3U,UniformCartesian<3U>,Vert >::iterator pf1, pf2, pf3;
        Index I = V1.getIndex(0);
        Index J = V1.getIndex(1);
        Index K = V1.getIndex(2);
        Vscalar1[I][J][K]= I;
        Vscalar2[I][J][K]= J;
        Vscalar3[I][J][K]= K;
        V1.Uncompress();
        for (pv=V1.begin(), pf1=Vscalar1.begin(),
                 pf2=Vscalar2.begin(), pf3=Vscalar3.begin() ;
             pv!=V1.end() ; ++pv, ++pf1, ++pf2, ++pf3) {
            (*pv)[0] = *pf1;
            (*pv)[1] = *pf2;
            (*pv)[2] = *pf3;
        }
    }
    //----------------------------------------------------------------------
    /*
    void initScalarVert(Field<double,1U,UniformCartesian<1U>,Vert>& V1) {
        Index I = V1.getIndex(0);
        V1[I]= I ;
    }
    //----------------------------------------------------------------------
    void initScalarVert(Field<double,2U,UniformCartesian<2U>,Vert>& V1) {
        Index I = V1.getIndex(0);
        Index J = V1.getIndex(1);
        V1[I][J]= I ;
    }
    */
    //----------------------------------------------------------------------
    void initScalarVert(Field<double,3U,UniformCartesian<3U>,Vert>& V1) {
        Index I = V1.getIndex(0);
        Index J = V1.getIndex(1);
        Index K = V1.getIndex(2);
        V1[I][J][K]= I ;
    }
}

TEST(Field,CellToVertex2D)
{
    OpalTestUtilities::SilenceTest silencer;

    const unsigned Dim = 2;
    int size = 8;
    NDIndex<Dim> cellDomain, vertDomain;
    for (unsigned int i = 0 ; i < Dim ; i++ ) {
        cellDomain[i] = Index(size);
        vertDomain[i] = Index(size+1);
    }

    // Some boundary condition objects:
    BConds<double,Dim,UniformCartesian<Dim>,Vert> bcScalarVert;
    BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalarCell;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcScalarVert[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Vert>(d);
        bcScalarCell[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
    }
    BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> bcVectorVert;
    BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVectorCell;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcVectorCell[d] =
            new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>(d);
        bcVectorVert[d] =
            new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>(d);
    }

    // perform the cell->vertex divergenace operation
    FieldLayout<Dim> layoutCell(cellDomain), layoutVert(vertDomain);
    Field<double,Dim,UniformCartesian<Dim>,Vert>
        ScaV1(layoutVert, bcScalarVert, GuardCellSizes<Dim>(1));
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>
        VecC1(layoutCell, bcVectorCell, GuardCellSizes<Dim>(1));

    initVectorCell(VecC1);

    // VecC1.writeb("vtest1");
    VecC1.fillGuardCells();
    // VecC1.writeb("vtest2");

    Inform* fdip =
        new Inform(NULL,"text.test.TestFieldCellToVertex2D",Inform::OVERWRITE,0);
    Inform& fdi = *fdip;
    setInform(fdi);
    //    setFormat(3,1); // 3 elements per line, 1 digit past decimal
    std::cout << "+++++++++VecC1 BEFORE begin+++++++++" << std::endl;
    fp2(VecC1);
    std::cout << "+++++++++VecC1 BEFORE end+++++++++++" << std::endl;
    std::cout << "sum VecC1 = " << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],224,roundOffError);

    //won't work w/new mesh classes (yet) --tjw assign(ScaV1 , Div(VecC1));
    ScaV1 = Div(VecC1,ScaV1);
    std::cout << "+++++++++ScaV1 = Div(VecC1,ScaV1) begin+++++++++" << std::endl;
    fp2(ScaV1);
    std::cout << "+++++++++ScaV1 = Div(VecC1,ScaV1) end+++++++++++" << std::endl;
    std::cout << "sum ScaV1 = " << sum(ScaV1) << std::endl;
    EXPECT_NEAR(sum(ScaV1),0,roundOffError);

    // ScaV1.writeb("ScaV1");

    // now perform the cell->vertex gradient operation
    Field<double,Dim,UniformCartesian<Dim>,Cell>
        ScaC1(layoutCell, bcScalarCell, GuardCellSizes<Dim>(1));
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>
        VecV1(layoutVert, bcVectorVert, GuardCellSizes<Dim>(1));

    initScalarCell(ScaC1);

    // ScaC1.writeb("stest1");
    ScaC1.fillGuardCells();
    // ScaC1.writeb("stest2");

    std::cout << "+++++++++ScaC1 BEFORE begin+++++++++" << std::endl;
    fp2(ScaC1);
    std::cout << "+++++++++ScaC1 BEFORE end+++++++++++" << std::endl;
    std::cout << "sum ScaC1 = " << sum(ScaC1) << std::endl;
    EXPECT_NEAR(sum(ScaC1),224,roundOffError);

    VecV1 = Grad(ScaC1,VecV1);
    std::cout << "+++++++++VecV1 = Grad(ScaC1,VecV1) begin+++++++++" << std::endl;
    fp2(VecV1);
    std::cout << "+++++++++VecV1 = Grad(ScaC1,VecV1) end+++++++++++" << std::endl;
    std::cout << "sum VecV1 = " << sum(VecV1) << std::endl;
    EXPECT_NEAR(sum(VecV1)[0],0,roundOffError);
    // VecV1.writeb("VecV1");
}

TEST(Field,CellToVertex3D)
{
    OpalTestUtilities::SilenceTest silencer;

    const unsigned Dim = 3;
    int size = 8;
    NDIndex<Dim> cellDomain;
    for (unsigned int i = 0 ; i < Dim ; i++ ) {
        cellDomain[i] = Index(size);
    }

    // perform the cell->vertex divergence operation
    FieldLayout<Dim> layoutCell(cellDomain);
    BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalar;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcScalar[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
    }
    Field<double,Dim,UniformCartesian<Dim>,Cell>
        ScaC1(layoutCell, bcScalar, GuardCellSizes<Dim>(1));
    BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVector;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcVector[d] =
            new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>(d);
    }
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>
        VecC1(layoutCell, bcVector, GuardCellSizes<Dim>(1));

    initVectorCell(VecC1);
    std::cout << VecC1 << std::endl;
    std::cout << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],1792,roundOffError);
    // VecC1.writeb("vtest1");

    VecC1.fillGuardCells();
    std::cout << VecC1 << std::endl;
    std::cout << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],1792,roundOffError);
    // VecC1.writeb("vtest2");

    ScaC1 = Div(VecC1,ScaC1);
    ScaC1.fillGuardCells();
    std::cout << ScaC1 << std::endl;
    std::cout << sum(ScaC1) << std::endl;
    EXPECT_NEAR(sum(ScaC1),672,roundOffError);
    // ScaC1.writeb("ScaC1");

    // now perform the cell->vertex gradient operation
    initScalarCell(ScaC1);
    std::cout << ScaC1 << std::endl;
    std::cout << sum(ScaC1) << std::endl;
    EXPECT_NEAR(sum(ScaC1),1792,roundOffError);
    // ScaC1.writeb("stest1");

    ScaC1.fillGuardCells();
    std::cout << ScaC1 << std::endl;
    std::cout << sum(ScaC1) << std::endl;
    EXPECT_NEAR(sum(ScaC1),1792,roundOffError);
    // ScaC1.writeb("stest2");

    VecC1 = Grad(ScaC1,VecC1);
    std::cout << VecC1 << std::endl;
    std::cout << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],512,roundOffError);
    // VecC1.writeb("VecC1");
}

TEST(Field, CelllToTenzorCell)
{
    OpalTestUtilities::SilenceTest silencer;
    const unsigned Dim = 2;
    int size = 4;

    NDIndex<Dim> cellDomain, vertDomain;
    for (unsigned int i = 0 ; i < Dim ; i++ ) {
        cellDomain[i] = Index(size);
        vertDomain[i] = Index(size+1);
    }

    FieldLayout<Dim> layoutCell(cellDomain), layoutVert(vertDomain);

    // perform the cell->tenzor divergence operation
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>
        VecV1(layoutCell);
    Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>
        TzrC1(layoutCell, GuardCellSizes<Dim>(1));

    initTenzorCell(TzrC1);
    std::cout << TzrC1 << std::endl;
    std::cout << sum(TzrC1) << std::endl;
    EXPECT_NEAR(sum(TzrC1)[0],24,roundOffError);
    // TzrC1.writeb("ttest1");

    TzrC1.fillGuardCells();
    std::cout << TzrC1 << std::endl;
    std::cout << sum(TzrC1) << std::endl;
    EXPECT_NEAR(sum(TzrC1)[0],24,roundOffError);
    // TzrC1.writeb("ttest2");

    Div(TzrC1,VecV1);
    std::cout << VecV1 << std::endl;
    std::cout << sum(VecV1) << std::endl;
    EXPECT_NEAR(sum(VecV1)[0],10.5,roundOffError);
    EXPECT_NEAR(sum(VecV1)[1],4.5,roundOffError);
    // VecV1.writeb("VecV1");

    // now perform the cell->tenzor gradient operation
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>
        VecC1(layoutCell, GuardCellSizes<Dim>(1));
    Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>
        TzrV1(layoutVert);

    initVectorCell(VecC1);
    std::cout << VecC1 << std::endl;
    std::cout << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],24,roundOffError);
    // VecC1.writeb("vtest1");

    VecC1.fillGuardCells();
    std::cout << VecC1 << std::endl;
    std::cout << sum(VecC1) << std::endl;
    EXPECT_NEAR(sum(VecC1)[0],24,roundOffError);
    // VecC1.writeb("vtest2");

    Grad(VecC1,TzrV1);
    std::cout << TzrV1 << std::endl;
    std::cout << sum(TzrV1) << std::endl;
    EXPECT_NEAR(sum(TzrV1)[0],0,roundOffError);
    // TzrV1.writeb("TrzV1");
}

TEST(Field,VertexToCell)
{
    OpalTestUtilities::SilenceTest silencer;

    const unsigned Dim = 3;
    int size = 8;

    NDIndex<Dim> cellDomain, vertDomain;
    for (unsigned int i = 0 ; i < Dim ; i++ ) {
        cellDomain[i] = Index(size);
        vertDomain[i] = Index(size+1);
    }

    // Some boundary condition objects:
    BConds<double,Dim,UniformCartesian<Dim>,Vert> bcScalarVert;
    BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalarCell;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcScalarVert[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Vert>(d);
        bcScalarCell[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
    }
    BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> bcVectorVert;
    BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVectorCell;
    for (unsigned int d=0; d<2*Dim; d++) {
        bcVectorCell[d] =
            new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>(d);
        bcVectorVert[d] =
            new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>(d);
    }

    // perform the vertex->cell divergenace operation
    FieldLayout<Dim> layoutCell(cellDomain), layoutVert(vertDomain);
    Field<double,Dim,UniformCartesian<Dim>,Cell>
        ScaC1(layoutCell, bcScalarCell, GuardCellSizes<Dim>(1));
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>
        VecV1(layoutVert, bcVectorVert, GuardCellSizes<Dim>(1));

    initVectorVert(VecV1);
    std::cout << VecV1 << std::endl;
    std::cout << sum(VecV1) << std::endl;
    EXPECT_NEAR(sum(VecV1)[0],2916,roundOffError);
    // VecV1.writeb("vtest1");

    VecV1.fillGuardCells();
    std::cout << VecV1 << std::endl;
    std::cout << sum(VecV1) << std::endl;
    EXPECT_NEAR(sum(VecV1)[0],2916,roundOffError);
    // VecV1.writeb("vtest2");

    ScaC1 = Div(VecV1,ScaC1);
    std::cout << ScaC1 << std::endl;
    std::cout << sum(ScaC1) << std::endl;
    EXPECT_NEAR(sum(ScaC1),1536,roundOffError);
    // ScaC1.writeb("ScaC1");

    // now perform the vertex->cell gradient operation
    Field<double,Dim,UniformCartesian<Dim>,Vert>
        ScaV1(layoutVert, bcScalarVert, GuardCellSizes<Dim>(1));
    Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>
        VecC1(layoutCell, bcVectorCell, GuardCellSizes<Dim>(1));

    initScalarVert(ScaV1);
    std::cout << ScaV1 << std::endl;
    std::cout << sum(ScaV1) << std::endl;
    EXPECT_NEAR(sum(ScaV1),2916,roundOffError);
    //ScaV1.writeb("stest1");

    ScaV1.fillGuardCells();
    std::cout << ScaV1 << std::endl;
    std::cout << sum(ScaV1) << std::endl;
    EXPECT_NEAR(sum(ScaV1),2916,roundOffError);
    //ScaV1.writeb("stest2");

    VecC1 = Grad(ScaV1,VecC1);
    std::cout << VecV1 << std::endl;
    std::cout << sum(VecV1) << std::endl;
    EXPECT_NEAR(sum(VecV1)[0],2916,roundOffError);
    //VecC1.writeb("VecC1");
}