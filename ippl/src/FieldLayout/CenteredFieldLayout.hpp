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

#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/Centering.h"
#include "Meshes/CartesianCentering.h"
#include "Utility/PAssert.h"


//=============================================================================
// Helper global functions:
// The constructors call these specialized global functions as a workaround for
// lack of partial specialization:
//=============================================================================

//===========================Arbitrary mesh type=============================

//-----------------------------------------------------------------------------
// These specify only a total number of vnodes, allowing the constructor
// complete control on how to do the vnode partitioning of the index space:
// Constructor for arbitrary dimension with parallel/serial specifier array:

//------------------Cell centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Cell> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++)
    ndi[d] = Index(mesh.gridSizes[d] - 1);
  cfl.initialize(ndi, edt, vnodes);
}

//=============================================================================
// General ctor calls specializations of global function (workaround for lack
// of partial specialization:
//=============================================================================

//-----------------------------------------------------------------------------
// These specify only a total number of vnodes, allowing the constructor
// complete control on how to do the vnode partitioning of the index space:
// Constructor for arbitrary dimension with parallel/serial specifier array:

// Constructor for arbitrary dimension with parallel/serial specifier array:
// This one also works if nothing except mesh is specified:
template<unsigned Dim, class Mesh, class Centering>
CenteredFieldLayout<Dim,Mesh,Centering>::
CenteredFieldLayout(Mesh& mesh,
		    e_dim_tag *p,
		    int vnodes)
{

  PInsist(Dim<=Mesh::Dimension,
    "CenteredFieldLayout dimension cannot be greater than Mesh dimension!!");
  centeredInitialize(*this, mesh, p, vnodes);
}

//-----------------------------------------------------------------------------
// These specify both the total number of vnodes and the numbers of vnodes
// along each dimension for the partitioning of the index space. Obviously this
// restricts the number of vnodes to be a product of the numbers along each
// dimension (the constructor implementation checks this):

// Constructor for arbitrary dimension with parallel/serial specifier array:
template<unsigned Dim, class Mesh, class Centering>
CenteredFieldLayout<Dim,Mesh,Centering>::
CenteredFieldLayout(Mesh& mesh,
		    e_dim_tag *p,
		    unsigned* vnodesAlongDirection,
		    bool recurse,
		    int vnodes)
{

  PInsist(Dim<=Mesh::Dimension,
    "CenteredFieldLayout dimension cannot be greater than Mesh dimension!!");
  centeredInitialize(*this, mesh, p, vnodesAlongDirection, recurse, vnodes);
}

//-----------------------------------------------------------------------------
// A constructor for a completely user-specified partitioning of the
// mesh space.

template<unsigned Dim, class Mesh, class Centering>
CenteredFieldLayout<Dim,Mesh,Centering>::
CenteredFieldLayout(Mesh& mesh,
		    const NDIndex<Dim> *dombegin,
		    const NDIndex<Dim> *domend,
		    const int *nbegin,
		    const int *nend)
{

  centeredInitialize(*this, mesh, dombegin, domend, nbegin, nend);
}
