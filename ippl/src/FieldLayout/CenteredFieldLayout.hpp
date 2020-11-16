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
#include "Meshes/CartesianCentering.h"
#include "Meshes/Centering.h"
#include "Utility/PAssert.h"

//-----------------------------------------------------------------------------
// These specify only a total number of vnodes, allowing the constructor
// complete control on how to do the vnode partitioning of the index space:
// Constructor for arbitrary dimension with parallel/serial specifier array:

//------------------Cell centering---------------------------------------------
template <unsigned Dim, class Mesh>
inline void centeredInitialize(
    CenteredFieldLayout<Dim, Mesh, Cell>& cfl, const Mesh& mesh, e_dim_tag* edt, int vnodes) {
    NDIndex<Dim> ndi;
    for (unsigned int d = 0; d < Dim; d++)
        ndi[d] = Index(mesh.gridSizes[d] - 1);
    cfl.initialize(ndi, edt, vnodes);
}

//------------------Vert centering---------------------------------------------
template <unsigned Dim, class Mesh>
inline void centeredInitialize(
    CenteredFieldLayout<Dim, Mesh, Vert>& cfl, const Mesh& mesh, e_dim_tag* edt, int vnodes) {
    NDIndex<Dim> ndi;
    for (unsigned int d = 0; d < Dim; d++)
        ndi[d] = Index(mesh.gridSizes[d]);
    cfl.initialize(ndi, edt, vnodes);
}

//-----------------------------------------------------------------------------
// These specify only a total number of vnodes, allowing the constructor
// complete control on how to do the vnode partitioning of the index space:
// Constructor for arbitrary dimension with parallel/serial specifier array:

// Constructor for arbitrary dimension with parallel/serial specifier array:
// This one also works if nothing except mesh is specified:
template <unsigned Dim, class Mesh, class Centering>
CenteredFieldLayout<Dim, Mesh, Centering>::CenteredFieldLayout(
    Mesh& mesh, e_dim_tag* p, int vnodes) {
    PInsist(
        Dim <= Mesh::Dimension,
        "CenteredFieldLayout dimension cannot be greater than Mesh dimension!!");
    centeredInitialize(*this, mesh, p, vnodes);
}
