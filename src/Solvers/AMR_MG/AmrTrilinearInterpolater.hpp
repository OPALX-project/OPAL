//
// Class AmrTrilinearInterpolater
//   Trilinear interpolation of data on coarse cells to fine cells.
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
AmrTrilinearInterpolater<Level>::AmrTrilinearInterpolater()
    : AmrInterpolater<Level>(2 << (AMREX_SPACEDIM - 1))
{ }


template <class Level>
void AmrTrilinearInterpolater<Level>::stencil(
    const AmrIntVect_t& iv,
    const basefab_t& fab,
    umap_t& map,
    const scalar_t& scale,
    Level* mglevel)
{
    /* lower left coarse cell (i, j, k)
     * floor( i - 0.5 ) / rr[0]
     * floor( j - 0.5 ) / rr[1]
     * floor( k - 0.5 ) / rr[2]
     */
    AmrIntVect_t civ;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            
        scalar_t tmp = iv[d] - 0.5;
        if ( std::signbit(tmp) )
            civ[d] = std::floor(tmp);
        else
            civ[d] = tmp;
    }
    
    civ.coarsen(mglevel->refinement());
        
    // ref ratio 2 only
    scalar_t dx = 0.5 * ( iv[0] - civ[0] * 2 ) - 0.25;
    scalar_t dy = 0.5 * ( iv[1] - civ[1] * 2 ) - 0.25;
#if AMREX_SPACEDIM == 3
    scalar_t dz = 0.5 * ( iv[2] - civ[2] * 2 ) - 0.25;
#endif
        
    scalar_t xdiff = 1.0 - dx;
    scalar_t ydiff = 1.0 - dy;
#if AMREX_SPACEDIM == 3
    scalar_t zdiff = 1.0 - dz;
#endif
    // (i, j, k)
    go_t crse_gidx = mglevel->serialize(civ);
    scalar_t value = AMREX_D_TERM(xdiff, * ydiff, * zdiff) * scale;
    
    if ( !mglevel->applyBoundary(civ, fab, map, value) )
        map[crse_gidx] += value;
    
    // (i+1, j, k)
    AmrIntVect_t tmp(D_DECL(civ[0]+1, civ[1], civ[2]));
    value = AMREX_D_TERM(dx, * ydiff, * zdiff) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
    
    // (i, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]));
    value = AMREX_D_TERM(xdiff, * dy, * zdiff) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
    
    // (i+1, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]));
    value = AMREX_D_TERM(dx, * dy, * zdiff) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
        
#if AMREX_SPACEDIM == 3
    // (i, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1], civ[2]+1));
    value = AMREX_D_TERM(xdiff, * ydiff, * dz) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
    
    // (i+1, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1], civ[2]+1));
    value = AMREX_D_TERM(dx, * ydiff, * dz) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
    
    // (i, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(xdiff, * dy, * dz) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
    
    // (i+1, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(dx, * dy, * dz) * scale;
    if ( !mglevel->applyBoundary(tmp, map, value) )
        map[mglevel->serialize(tmp)] += value;
#endif
}


template <class Level>
void AmrTrilinearInterpolater<Level>::coarse(
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
void AmrTrilinearInterpolater<Level>::fine(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& fab,
    Level* mglevel)
{
    /*
     * The AmrTrilinearInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, fab, map, scale, mglevel);
}
