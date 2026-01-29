//
// Class AmrBoundary
//   Interface for mesh boundary conditions.
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
#ifndef AMR_BOUNDARY_H
#define AMR_BOUNDARY_H

#include <vector>

#include "Amr/AmrDefs.h"

template <class Level>
class AmrBoundary {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef typename Level::basefab_t   basefab_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    /*!
     * @param nPoints used in stencil for applying the boundary
     */
    AmrBoundary(go_t nPoints) : nPoints_m(nPoints) { };
    
    
    /*!
     * Check if we are on the physical boundary (all directions)
     * @param iv cell to check
     * @param nr is the number of grid points
     */
    bool isBoundary(const AmrIntVect_t& iv, const go_t* nr) const {
        return AMREX_D_TERM(   isBoundary(iv, 0, nr),
                            || isBoundary(iv, 1, nr),
                            || isBoundary(iv, 2, nr));
    }

    /*!
     * Check if we are on the physical boundary (certain direction)
     * @param iv cell to check
     * @param nr is the number of grid points
     */
    bool isBoundary(const AmrIntVect_t& iv,
                    const lo_t& dir,
                    const go_t* nr) const {
        return ( iv[dir] < 0 || iv[dir] >= nr[0] );
    }
    
    /*!
     * Apply boundary in a certain direction.
     * @param iv is the cell where we want to have the boundary value
     * @param dir direction of physical / mesh boundary
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     * @param value of matrix entry that is supposed for index
     * @param nr is the number of grid points
     */
    virtual void apply(const AmrIntVect_t& iv,
                       const lo_t& dir,
                       umap_t& map,
                       const scalar_t& value,
                       Level* mglevel,
                       const go_t* nr) = 0;
    
    /*!
     * @returns the number of stencil points required
     */
    const go_t& getNumberOfPoints() const {
        return nPoints_m;
    }
    
private:
    const go_t nPoints_m;    ///< Number of points used for boundary
};

#endif
