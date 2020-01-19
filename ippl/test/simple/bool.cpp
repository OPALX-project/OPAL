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
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "Message/Communicate.h"

#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>

// forward declarations
ofstream *logfile;


int main(int argc, char *argv[] )
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  char name[100];
  sprintf(name,"out.%d",Ippl::Comm->myNode());
  logfile = new ofstream(name);
  *logfile << "Starting to run!" << endl;
  Ippl::Comm->barrier();

  Index I(10);
  const unsigned Dim=1;

  double amin;
  FieldLayout<Dim>  layout(I);
  Field<double,Dim,UniformCartesian<Dim>,
    UniformCartesian<Dim>::DefaultCentering> A(layout);
  Field<bool,Dim,UniformCartesian<Dim>,
    UniformCartesian<Dim>::DefaultCentering>   B(layout);
  Field<double,Dim,UniformCartesian<Dim>,
    UniformCartesian<Dim>::DefaultCentering> C(layout);

  A[I] = I;
  amin = 5.0;
  C = amin;
  B = lt( A, amin ) ;
  A = where( B, C, A ) ;
  A = where( lt(A,amin), A, C ) ;

  return 0;
}
