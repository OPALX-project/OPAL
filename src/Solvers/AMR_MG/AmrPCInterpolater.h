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
#ifndef AMR_PIECEWISE_CONST_INTERPOLATER_H
#define AMR_PIECEWISE_CONST_INTERPOLATER_H

#include "AmrInterpolater.h"

template <class Level>
class AmrPCInterpolater : public AmrInterpolater<Level>
{
public:
    typedef typename Level::go_t        go_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::basefab_t   basefab_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    AmrPCInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 const basefab_t& fab,
                 umap_t& map,
                 const scalar_t& scale,
                 Level* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const basefab_t& rfab,
                const AmrIntVect_t& riv,
                Level* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
              const scalar_t& scale,
              lo_t dir, lo_t shift, const basefab_t& fab,
              Level* mglevel);
    
};

#include "AmrPCInterpolater.hpp"

#endif
