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

//
// Print out an array including the guard cells.
//

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

//-----------------------------------------------------------------------------
// class LogicalY
//
// A simple demonstration of setting spacially dependent boundary
// conditions.  It inherits from PatchBC<T,D,M,C> and defines an apply
// member function which applies the boundary condition to a given
// domain using a Field iterator.
//-----------------------------------------------------------------------------

template <class T, unsigned D, class M, class C>
class LogicalY : public PatchBC<T,D,M,C>
{
public:

  //
  // Construct the boundary condition using an integer for the face.
  // The format of this integer is the same as for any other boundary
  // condition.  It is stored in the base class, and you can get it by
  // calling the member function getFace();
  //
  LogicalY(int f) : PatchBC<T,D,M,C>(f) {}

  //
  // Definition of a simple boundary condition.  In this case it just
  // sets the value in the boundary to the logical y coordinate.
  //
  void applyPatch(typename Field<T,D,M,C>::iterator p , const NDIndex<D>& x)
    {
      //
      // Get the size of the domain so we can loop over it.
      //
      int nx = x[0].length();
      int ny = x[1].length();

      //
      // Get the offset of the domain so we can use the global coordinate.
      //
      int ox = x[0].first();
      int oy = x[1].first();

      //
      // Loop over the domain setting each value.  Note that we use
      // the *local coordinate* as arguments to offset.  That is
      // because the iterator starts out pointed to the lower left
      // corner of the boundary domain.  On the right hand side we add
      // the offset back in so that we're setting the global
      // coordinate.
      //
      for (int x=0; x<nx; ++x)
	for (int y=0; y<ny; ++y)
	  p.offset(x,y) = y+oy;
    }

  //
  // The interface of boundary condition objects requires a clone
  // member function in the leaf classes.
  //
  BCondBase<T,D,M,C>* clone() const
    {
      return new LogicalY<T,D,M,C>( *this );
    }
};

//-----------------------------------------------------------------------------


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

  // Set initial boundary conditions.
  typedef LogicalY<double,Dim,M,C> TP;
  BConds<double,Dim,M,C> bc;
  bc[0] = new TP(0);
  bc[1] = new TP(1);
  bc[2] = new TP(2);
  bc[3] = new TP(3);

  // An array for testing.
  Field<double,Dim,M,C> A(layout,GuardCellSizes<Dim>(2),bc);

  // Fill A with zeros.
  A = 0.0;
  cout << "A:\n";
  printout(A);

  return 0;
}

/***************************************************************************
 * $RCSfile: bc3.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: bc3.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
