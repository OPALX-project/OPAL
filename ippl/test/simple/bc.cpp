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
#include "Field/GuardCellSizes.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"

#include <iostream>

// define helper print function
void printout(Field<double,2>& A)
{
  GuardCellSizes<2> gc = A.getGC();
  NDIndex<2> domain = A.getDomain();
  double *p = &*A.begin();
  int l0 = domain[0].length() + gc.left(0) + gc.right(0);
  int l1 = domain[1].length() + gc.left(1) + gc.right(1);
  p -= gc.left(0);
  p -= gc.left(1)*l0;

  for (int j=0; j<l1; ++j)
    {
      for (int i=0; i<l0; ++i, ++p)
	cout << *p << " ";
      cout << endl;
    }
  cout << endl;
  return;
}


int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=2;
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
  testmsg << "A = " << A << endl;
  B[I][J] = A[I-2][J-2];
  testmsg << "B = " << B << endl;
  B[I][J] = A[I+2][J+2];
  testmsg << "B = " << B << endl;
  B[I][J] = A[I+2][J-2];
  testmsg << "B = " << B << endl;
  B[I][J] = A[I-2][J+2];
  testmsg << "B = " << B << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: bc.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: bc.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
