//
// Class AmrInterpolater
//   Abstract base class for all coarse to fine cell interpolaters.
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
#ifndef AMR_INTERPOLATER_H
#define AMR_INTERPOLATER_H

#include "Amr/AmrDefs.h"

template <class Level>
class AmrInterpolater {

public:
    typedef typename Level::go_t        go_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::basefab_t   basefab_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    /*!
     * @param nPoints is the number of interpolation points used
     */
    AmrInterpolater(lo_t nPoints) : nPoints_m(nPoints) { }
    
    /*!
     * Number of cell points used for interpolation.
     */
    const lo_t& getNumberOfPoints() const {
        return nPoints_m;
    }
    
    /*!
     * Get the stencil to interpolate a value from coarse to fine level
     * @param iv is the fine cell where we want to have the interpolated value
     * @param map with global matrix indices of coarse level cells and
     * matrix entries of coarse level cells (coefficients)
     * @param scale to apply to matrix values
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void stencil(const AmrIntVect_t& iv,
                         const basefab_t& fab,
                         umap_t& map,
                         const scalar_t& scale,
                         Level* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of coarse side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * matrix entries of coarse level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param riv is the fine cell at the interface
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary values at physical domain, e.g. Dirichlet, open BC
     */
    virtual void coarse(const AmrIntVect_t& iv,
                        umap_t& map,
                        const scalar_t& scale,
                        lo_t dir, lo_t shift, const basefab_t& rfab,
                        const AmrIntVect_t& riv,
                        Level* mglevel) = 0;
    
    /*!
     * Coarse-Fine-Interface
     * Get stencil of fine side
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param map with global matrix indices of fine level cells and
     * matrix entries of fine level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param ba contains all coarse cells that got refined
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    virtual void fine(const AmrIntVect_t& iv,
                      umap_t& map,
                      const scalar_t& scale,
                      lo_t dir, lo_t shift, const basefab_t& fab,
                      Level* mglevel)
    { };

    /*!
     * Coarse-Fine-Interface
     * Get stencil of fine side
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param map with global matrix indices of fine level cells and
     * matrix entries of fine level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    virtual void fine(const AmrIntVect_t& iv,
                      umap_t& map,
                      const scalar_t& scale,
                      lo_t dir, lo_t shift,
                      Level* mglevel)
    { };
    
protected:
    const lo_t nPoints_m;    ///< Number of points used for interpolation
};

#endif
