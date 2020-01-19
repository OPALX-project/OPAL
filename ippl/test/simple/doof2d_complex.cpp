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

int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);
  Inform msg("doof2d");

  int Parent = 0;
  int tag = 101;

  int vnodes, sizeX, sizeY, centerX, centerY, iterations;

  // read in data off of node 0
  if( Ippl::Comm->myNode() == Parent ) {
    cout << "input the number of vnodes" << endl;
    cin >> vnodes; 
    cout << "input sizeX" << endl;
    cin >> sizeX;
    cout << "input sizeY" << endl;
    cin >> sizeY;
    cout << "input centerX" << endl;
    cin >> centerX;
    cout << "input centerY" << endl;
    cin >> centerY;
    cout << "input the iterations" << endl;
    cin >> iterations;
    // now broadcast data to other nodes
    
    Message *mess = new Message();
    putMessage( *mess, vnodes );
    putMessage( *mess, sizeX );
    putMessage( *mess, sizeY );
    putMessage( *mess, centerX );
    putMessage( *mess, centerY );
    putMessage( *mess, iterations );
    Ippl::Comm->broadcast_all(mess, tag);
  }
  // now each node recieves the data
  
  Message *mess = Ippl::Comm->receive_block(Parent, tag);
  PAssert(mess);
  getMessage( *mess, vnodes );
  getMessage( *mess, sizeX );
  getMessage( *mess, sizeY );
  getMessage( *mess, centerX );
  getMessage( *mess, centerY );
  getMessage( *mess, iterations );
  delete mess;
   
  msg << "received a message on node " << Ippl::Comm->myNode();
  msg << " from node " << Parent << " with tag " << tag << endl;

  const unsigned Dim=2;
  Index I(sizeX);
  Index J(sizeY);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL, vnodes);
  Field<std::complex<double>,Dim> A(layout,GuardCellSizes<Dim>(1));
  Field<std::complex<double>,Dim> B(layout);

  A = 0.0;
  A[centerX][centerY] = std::complex<double>(1.0*iterations,0.0);

  double fact = 1.0/9.0;
  for(int iter = 0 ; iter < iterations ; iter++ ) 
  {
    assign(A[I][J], fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] + 
                          A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] + 
                          A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1]));
    B= 1.0*iter;
    cout << "iter = " << iter << " sum = " << sum(A) << endl;
  }
  cout << "enter to continue" << endl;
  int in;
  cin >> in;
  return 0;
}

/***************************************************************************
 * $RCSfile: doof2d_complex.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: doof2d_complex.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
