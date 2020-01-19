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

// TestPeriodicBC.cpp
// Cell centered periodic BC

// include files
#include "Ippl.h"
#include "Utility/FieldDebug.h"

#include <fstream>

typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  setFormat(3,1);

  const unsigned Dim = 3;
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

  msg << "++++++++BConds object cbc begin++++++++" << endl;
  msg << cbc;
  msg << "++++++++BConds object cbc end++++++++++" << endl;

  // Cell-centered Field's:
  msg << "layout: " << layout << endl;
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
  msg << "++++++++++cA+++++++++++" << endl ;
  fp3(cA);
  cB[I][J][K] = cA[I][J][K-2];
  msg << "++++++++++cB+++++++++++" << endl ;
  fp3(cB);

  dVf = Grad(cA,dVf);
  msg << "++++++++++Grad(cA)+++++" << endl ;
  fp3(dVf);

  dVf = Vector_t(0.0);

  dVf = Grad(cB,dVf);
  msg << "++++++++++Grad(cB)+++++" << endl ;
  fp3(dVf);



}
/***************************************************************************
 * $RCSfile: TestPeriodicBC.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: TestPeriodicBC.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
