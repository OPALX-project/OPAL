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

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(domain);
  Field<double,Dim> A(layout);
  Field<int,Dim> B(layout,GuardCellSizes<Dim>(1));

  Field<int,Dim>::iterator p;
  int i;
  A = -1.0 ;
  for (p=B.begin(), i=0; p!=B.end(); ++p, ++i) *p = i;

  testmsg << A << endl;
  testmsg << B << endl;

  A[I][J] = B[I-1][J+1] + B[I+1][J-1] ;
  testmsg << A << endl;

  B.fillGuardCells();
  A[I][J] = B[I-1][J+1] + B[I+1][J-1];
  testmsg << A << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: t.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: t.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
