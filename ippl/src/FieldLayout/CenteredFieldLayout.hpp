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

//------------------Vert centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Vert> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++)
    ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, edt, vnodes);
}

//------------------Edge centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Edge> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++)
    ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, edt, vnodes);
}

//------------------CartesianCentering centering-------------------------------
template<CenteringEnum* CE, unsigned Dim, class Mesh,
         unsigned NComponents>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,
		   CartesianCentering<CE,Dim,NComponents> > & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  // For componentwise layout of Field of multicomponent object, like
  // Field<Vektor...>, allow for maximal number of objects needed per
  // dimension (the number for the object component requiring the maximal
  // number):
  unsigned npts[Dim], nGridPts;
  unsigned int d, c;
  for (d=0; d<Dim; d++) {
    nGridPts = mesh.gridSizes[d];
    npts[d] = 0;
    for (c=0; c<NComponents; c++) {
      if (CE[c + d*NComponents] == CELL) {
        npts[d] = std::max(npts[d], (nGridPts - 1));
      } else {
        npts[d] = std::max(npts[d], nGridPts);
      }
    }
  }
  for (d=0; d<Dim; d++) ndi[d] = Index(npts[d]);
  cfl.initialize(ndi, edt, vnodes);
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// These specify both the total number of vnodes and the numbers of vnodes
// along each dimension for the partitioning of the index space. Obviously this
// restricts the number of vnodes to be a product of the numbers along each
// dimension (the constructor implementation checks this):

//------------------Cell centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Cell> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   unsigned* vnodesAlongDirection,
		   bool recurse,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d] - 1);
  cfl.initialize(ndi, edt, vnodesAlongDirection, recurse, vnodes);
}

//------------------Vert centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Vert> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   unsigned* vnodesAlongDirection,
		   bool recurse,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, edt, vnodesAlongDirection, recurse, vnodes);
}

//------------------Edge centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Edge> & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   unsigned* vnodesAlongDirection,
		   bool recurse,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, edt, vnodesAlongDirection, recurse, vnodes);
}

//------------------CartesianCentering centering-------------------------------
template<CenteringEnum* CE, unsigned Dim, class Mesh,
         unsigned NComponents>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,
		   CartesianCentering<CE,Dim,NComponents> > & cfl,
		   const Mesh& mesh,
		   e_dim_tag* edt,
		   unsigned* vnodesAlongDirection,
		   bool recurse,
		   int vnodes)
{
  NDIndex<Dim> ndi;
  // For componentwise layout of Field of multicomponent object, like
  // Field<Vektor...>, allow for maximal number of objects needed per
  // dimension (the number for the object component requiring the maximal
  // number):
  unsigned npts[Dim], nGridPts;
  unsigned int d, c;
  for (d=0; d<Dim; d++) {
    nGridPts = mesh.gridSizes[d];
    npts[d] = 0;
    for (c=0; c<NComponents; c++) {
      if (CE[c + d*NComponents] == CELL) {
	npts[d] = max(npts[d], (nGridPts - 1));
      } else {
	npts[d] = max(npts[d], nGridPts);
      }
    }
  }
  for (d=0; d<Dim; d++) ndi[d] = Index(npts[d]);
  cfl.initialize(ndi, edt, vnodesAlongDirection, recurse, vnodes);
}

//-----------------------------------------------------------------------------
// This initializer just specifies a completely user-provided partitioning.

//------------------Cell centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Cell> & cfl,
		   const Mesh& mesh,
		   const NDIndex<Dim> *dombegin,
		   const NDIndex<Dim> *domend,
		   const int *nbegin,
		   const int *nend)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d] - 1);
  cfl.initialize(ndi, dombegin, domend, nbegin, nend);
}

//------------------Vert centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Vert> & cfl,
		   const Mesh& mesh,
		   const NDIndex<Dim> *dombegin,
		   const NDIndex<Dim> *domend,
		   const int *nbegin,
		   const int *nend)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, dombegin, domend, nbegin, nend);
}

//------------------Edge centering---------------------------------------------
template<unsigned Dim, class Mesh>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,Edge> & cfl,
		   const Mesh& mesh,
		   const NDIndex<Dim> *dombegin,
		   const NDIndex<Dim> *domend,
		   const int *nbegin,
		   const int *nend)
{
  NDIndex<Dim> ndi;
  for (unsigned int d=0; d<Dim; d++) ndi[d] = Index(mesh.gridSizes[d]);
  cfl.initialize(ndi, dombegin, domend, nbegin, nend);
}

//------------------CartesianCentering centering-------------------------------
template<CenteringEnum* CE, unsigned Dim, class Mesh,
         unsigned NComponents>
inline void
centeredInitialize(CenteredFieldLayout<Dim,Mesh,
		   CartesianCentering<CE,Dim,NComponents> > & cfl,
		   const Mesh& mesh,
		   const NDIndex<Dim> *dombegin,
		   const NDIndex<Dim> *domend,
		   const int *nbegin,
		   const int *nend)
{
  // For componentwise layout of Field of multicomponent object, like
  // Field<Vektor...>, allow for maximal number of objects needed per
  // dimension (the number for the object component requiring the maximal
  // number):
  unsigned npts[Dim], nGridPts;
  unsigned int d, c;
  for (d=0; d<Dim; d++) {
    nGridPts = mesh.gridSizes[d];
    npts[d] = 0;
    for (c=0; c<NComponents; c++) {
      if (CE[c + d*NComponents] == CELL) {
	npts[d] = max(npts[d], (nGridPts - 1));
      } else {
	npts[d] = max(npts[d], nGridPts);
      }
    }
  }

  NDIndex<Dim> ndi;
  for (d=0; d<Dim; d++) ndi[d] = Index(npts[d]);
  cfl.initialize(ndi, dombegin, domend, nbegin, nend);
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
