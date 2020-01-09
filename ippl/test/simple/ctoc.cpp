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

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Ippl.h"

void initVectorCell(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1);
void initVectorCell(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1);
void initVectorCell(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1);
void initScalarCell(Field<double,1U,UniformCartesian<1U>,Cell>& C1);
void initScalarCell(Field<double,2U,UniformCartesian<2U>,Cell>& C1);
void initScalarCell(Field<double,3U,UniformCartesian<3U>,Cell>& C1);

const unsigned Dim = 3;
int size = 8;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  
  NDIndex<Dim> cellDomain;
  for (int i = 0 ; i < Dim ; i++ ) {
    cellDomain[i] = Index(size);
  }

  // perform the cell->vertex divergenace operation 
  FieldLayout<Dim> layoutCell(cellDomain);
  BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalar;
  for (int d=0; d<2*Dim; d++) {
    bcScalar[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
  }
  Field<double,Dim,UniformCartesian<Dim>,Cell>  
    ScaC1(layoutCell, bcScalar, GuardCellSizes<Dim>(1));
  BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVector;
  for (int d=0; d<2*Dim; d++) {
    bcVector[d] = 
      new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>(d);
  }
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> 
    VecC1(layoutCell, bcVector, GuardCellSizes<Dim>(1));
  
  initVectorCell(VecC1);

  VecC1.writeb("vtest1");
  VecC1.fillGuardCells();
  VecC1.writeb("vtest2");

  ScaC1 = Div(VecC1,ScaC1);
  ScaC1.fillGuardCells();
 
  ScaC1.writeb("ScaC1");

  // now perform the cell->vertex gradient operation 
  initScalarCell(ScaC1);

  ScaC1.writeb("stest1");
  ScaC1.fillGuardCells();
  ScaC1.writeb("stest2");

  VecC1 = Grad(ScaC1,VecC1);

  VecC1.writeb("VecC1");

  return 0;
}
//----------------------------------------------------------------------
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
void initScalarCell(Field<double,1U,UniformCartesian<1U>,Cell>& C1) {
  Index I = C1.getIndex(0);
  C1[I]= I ;
}
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

/***************************************************************************
 * $RCSfile: ctoc.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: ctoc.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
