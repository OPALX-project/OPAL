//
// Class CenteredFieldLayout
//
// Copyright (c) 2003 - 2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef CENTERED_FIELD_LAYOUT_H
#define CENTERED_FIELD_LAYOUT_H

#include "FieldLayout/FieldLayout.h"

template<unsigned Dim, class Mesh, class Centering>
class CenteredFieldLayout : public FieldLayout<Dim>
{
public:
  //---------------------------------------------------------------------------
  // Constructors from a mesh object only and parallel/serial specifiers.
  // If not doing this, user should be just using simple FieldLayout object, 
  // though no harm would be done in constructiong a CenteredFieldLayout with
  // Index/NDIndex arguments via the inherited constructors from FieldLayout.
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  // These specify only a total number of vnodes, allowing the constructor
  // complete control on how to do the vnode partitioning of the index space:
  // Constructor for arbitrary dimension with parallel/serial specifier array:

  // This one also works if nothing except mesh is specified:
  CenteredFieldLayout(Mesh& mesh, 
		      e_dim_tag *p=0, 
		      int vnodes=-1);

  // Special constructor which uses a existing partition
  // particular from expde


  //---------------------------------------------------------------------------
  // These specify both the total number of vnodes and the numbers of vnodes
  // along each dimension for the partitioning of the index space. Obviously
  // this restricts the number of vnodes to be a product of the numbers along
  // each dimension (the constructor implementation checks this):

  // Constructor for arbitrary dimension with parallel/serial specifier array:
  CenteredFieldLayout(Mesh& mesh, e_dim_tag *p, 
		      unsigned* vnodesAlongDirection, 
		      bool recurse=false,
		      int vnodes=-1);
  //---------------------------------------------------------------------------
  // A constructor a a completely user-specified partitioning of the
  // mesh space.

  CenteredFieldLayout(Mesh& mesh,
		      const NDIndex<Dim> *dombegin,
		      const NDIndex<Dim> *domend,
		      const int *nbegin,
		      const int *nend);
};

#include "FieldLayout/CenteredFieldLayout.hpp"

#endif // CENTERED_FIELD_LAYOUT_H
