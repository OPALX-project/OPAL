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
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/FieldSpec.h"
#include "Field/GuardCellSizes.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"

#include <iostream>

// set dimensionality
const unsigned Dim=2;

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


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

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
  testmsg << A << endl;
  A.writeb("atest");
  assign(B[I][J], A[I-1][J-1]);
  testmsg << B << endl;
  B.writeb("btest1");
  assign(B[I][J], A[I+1][J+1]);
  testmsg << B << endl;
  B.writeb("btest2");

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
  testmsg << A2 << endl;
  A.writeb("a2test");
  assign(B[I][J], A[I-1][J-1]);
  testmsg << B << endl;
  B.writeb("b2test1");
  assign(B[I][J], A[I+1][J+1]);
  testmsg << B << endl;
  B.writeb("b2test2");

  return 0;
}

/***************************************************************************
 * $RCSfile: bc2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: bc2.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
