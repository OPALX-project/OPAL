//
// Class AmrPCInterpolater
//   Piecewise constant interpolation of data on coarse cells to fine cells.
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
template <class Level>
AmrPCInterpolater<Level>::AmrPCInterpolater()
    : AmrInterpolater<Level>(1)
{ }


template <class Level>
void AmrPCInterpolater<Level>::stencil(
    const AmrIntVect_t& iv,
    const basefab_t& fab,
    umap_t& map,
    const scalar_t& scale,
    Level* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    if ( !mglevel->applyBoundary(civ, fab, map, scale) )
        map[mglevel->serialize(civ)] += scale;
}


template <class Level>
void AmrPCInterpolater<Level>::coarse(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    Level* mglevel)
{
    // do nothing
}


template <class Level>
void AmrPCInterpolater<Level>::fine(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& fab,
    Level* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, fab, map, scale, mglevel);
}
