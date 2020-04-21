//
// Class AmrDirichletBoundary
//   Dirichlet boundary is on faces of physical domain the boundary
//   value would be at different locations depending on the level.
//
// Copyright (c) 2017 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef AMR_DIRICHLET_BOUNDARY_H
#define AMR_DIRICHLET_BOUNDARY_H

#include "AmrBoundary.h"

template <class Level>
class AmrDirichletBoundary : public AmrBoundary<Level> {
    
public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    AmrDirichletBoundary() : AmrBoundary<Level>(1) { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
};


template <class Level>
void AmrDirichletBoundary<Level>::apply(const AmrIntVect_t& iv,
                                        const lo_t& dir,
                                        umap_t& map,
                                        const scalar_t& value,
                                        Level* mglevel,
                                        const go_t* nr)
{
    // find interior neighbour cell
    AmrIntVect_t niv = iv;
    niv[dir] = (iv[dir] == -1) ? iv[dir] + 1 : iv[dir] - 1;
    map[mglevel->serialize(niv)] -= value;
}




#endif
