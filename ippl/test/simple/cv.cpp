// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include <iostream>
#include "Ippl.h"

void initVectorCell(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1);
void initVectorCell(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1);
void initVectorCell(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1);
void initScalarCell(Field<double,1U,UniformCartesian<1U>,Cell>& C1);
void initScalarCell(Field<double,2U,UniformCartesian<2U>,Cell>& C1);
void initScalarCell(Field<double,3U,UniformCartesian<3U>,Cell>& C1);

const unsigned Dim = 2;
int size = 8;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  
  NDIndex<Dim> cellDomain, vertDomain;
  for (int i = 0 ; i < Dim ; i++ ) {
    cellDomain[i] = Index(size);
    vertDomain[i] = Index(size+1);
  }

  // Some boundary condition objects:
  BConds<double,Dim,UniformCartesian<Dim>,Vert> bcScalarVert;
  BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalarCell;
  for (int d=0; d<2*Dim; d++) {
    bcScalarVert[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Vert>(d);
    bcScalarCell[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
  }
  BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> bcVectorVert;
  BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVectorCell;
  for (int d=0; d<2*Dim; d++) {
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

  VecC1.writeb("vtest1");
  VecC1.fillGuardCells();
  VecC1.writeb("vtest2");

  setFormat(3,1); // 3 elements per line, 1 digit past decimal
//   cout << "+++++++++VecC1 BEFORE begin+++++++++" << endl;
//   fp2(VecC1);
//   cout << "+++++++++VecC1 BEFORE end+++++++++++" << endl;
  //won't work w/new mesh classes (yet) --tjw assign(ScaV1 , Div(VecC1));
  ScaV1 = Div(VecC1,ScaV1);
//   cout << "+++++++++ScaV1 = Div(VecC1,ScaV1) begin+++++++++" << endl;
//   fp2(ScaV1);
//   cout << "+++++++++ScaV1 = Div(VecC1,ScaV1) end+++++++++++" << endl;
 
  ScaV1.writeb("ScaV1");

  // now perform the cell->vertex gradient operation 
  Field<double,Dim,UniformCartesian<Dim>,Cell> 
    ScaC1(layoutCell, bcScalarCell, GuardCellSizes<Dim>(1));
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> 
    VecV1(layoutVert, bcVectorVert, GuardCellSizes<Dim>(1));

  initScalarCell(ScaC1);

  ScaC1.writeb("stest1");
  ScaC1.fillGuardCells();
  ScaC1.writeb("stest2");

//   cout << "+++++++++ScaC1 BEFORE begin+++++++++" << endl;
//   fp2(ScaC1);
//   cout << "+++++++++ScaC1 BEFORE end+++++++++++" << endl;
  VecV1 = Grad(ScaC1,VecV1);
//   cout << "+++++++++VecV1 = Grad(ScaC1,VecV1) begin+++++++++" << endl;
//   fp2(VecV1);
//   cout << "+++++++++VecV1 = Grad(ScaC1,VecV1) end+++++++++++" << endl;

  VecV1.writeb("VecV1");

  return 0;
}
//----------------------------------------------------------------------
void 
initVectorCell(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1) {
  //1.0 is default spacing  double X[1] = { 1.0 } ;
  //  UniformCartesian<1U>::setMesh(X);
  //1.0 is default spacing  C1.get_mesh().set_meshSpacing(X);
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
//----------------------------------------------------------------------
void 
initVectorCell(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1) {
  //1.0 is default spacing  double X[2] = { 1.0, 1.0 } ;
  //  UniformCartesian<2U>::setMesh(X);
  //1.0 is default spacing  C1.get_mesh().set_meshSpacing(X);
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
  //1.0 is default spacing  double X[3] = { 1.0, 1.0, 1.0 } ;
  //  UniformCartesian<3U>::setMesh(X);
  //1.0 is default spacing  C1.get_mesh().set_meshSpacing(X);
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
void initScalarCell(Field<double,1U,UniformCartesian<1U>,Cell>& V1) {
  Index I = V1.getIndex(0);
  V1[I]= I ;
}
//----------------------------------------------------------------------
void initScalarCell(Field<double,2U,UniformCartesian<2U>,Cell>& V1) {
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  V1[I][J]= I ;
}
//----------------------------------------------------------------------
void initScalarCell(Field<double,3U,UniformCartesian<3U>,Cell>& V1) {
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  Index K = V1.getIndex(2);
  V1[I][J][K]= I ;
}
//----------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: cv.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: cv.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
