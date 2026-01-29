//
// Class AmrOpenBoundary
//   Open boundary condition
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

#ifndef AMR_OPEN_BOUNDARY_H
#define AMR_OPEN_BOUNDARY_H

#include "AmrBoundary.h"

#include "Utilities/OpalException.h"

template <class Level>
class AmrOpenBoundary : public AmrBoundary<Level> {

public:
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::go_t        go_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
public:
    
    AmrOpenBoundary()
        : AmrBoundary<Level>(4)
        , order_m(ABC::Robin)
        , dist_m(1.7)
    { }
    
    void apply(const AmrIntVect_t& iv,
               const lo_t& dir,
               umap_t& map,
               const scalar_t& value,
               Level* mglevel,
               const go_t* nr);
    
    enum ABC {
        Zeroth,
        First,
        Second,
        Third,
        Robin
    };
    
private:
    int order_m;
    double dist_m;
    
    scalar_t coordinate_m(const AmrIntVect_t& iv,
                          const lo_t& dir,
                          Level* mglevel,
                          const go_t* nr);
    
    /*!
     * Robin boundary condition
     */
    void robin_m(const AmrIntVect_t& iv,
                 const lo_t& dir,
                 umap_t& map,
                 const scalar_t& value,
                 Level* mglevel,
                 const go_t* nr);
    
    /*!
     * Asymptotic boundary condition 0th order (ABC0)
     */
    void abc0_m(const AmrIntVect_t& iv,
                const lo_t& dir,
                umap_t& map,
                const scalar_t& value,
                Level* mglevel,
                const go_t* nr);
    
    
    /*!
     * Asymptotic boundary condition 1st order (ABC1)
     */
    void abc1_m(const AmrIntVect_t& iv,
                const lo_t& dir,
                umap_t& map,
                const scalar_t& value,
                Level* mglevel,
                const go_t* nr);
    
    /*!
     * Asymptotic boundary condition 2nd order (ABC2)
     */
    void abc2_m(const AmrIntVect_t& iv,
                const lo_t& dir,
                umap_t& map,
                const scalar_t& value,
                Level* mglevel,
                const go_t* nr);
    
    
    /*!
     * Asymptotic boundary condition 3rd order (ABC3)
     */
    void abc3_m(const AmrIntVect_t& iv,
                const lo_t& dir,
                umap_t& map,
                const scalar_t& value,
                Level* mglevel,
                const go_t* nr);
    
};


template <class Level>
void AmrOpenBoundary<Level>::apply(const AmrIntVect_t& iv,
                                   const lo_t& dir,
                                   umap_t& map,
                                   const scalar_t& value,
                                   Level* mglevel,
                                   const go_t* nr)
{
    switch ( order_m ) {
        case ABC::Third:
        {
            this->abc3_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Second:
        {
            this->abc2_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::First:
        {
            this->abc1_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Zeroth:
        {
            this->abc0_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        case ABC::Robin:
        {
            this->robin_m(iv, dir, map, value, mglevel, nr);
            break;
        }
        default:
        {
            throw OpalException("AmrOpenBoundary<Level>::apply()",
                                "Order not supported.");
        }
    }
}


template <class Level>
void AmrOpenBoundary<Level>::robin_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     umap_t& map,
                                     const scalar_t& value,
                                     Level* mglevel,
                                     const go_t* nr)
{
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        n2iv[dir] = 1;
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    // correct cell size
    scalar_t h = mglevel->cellSize(dir);
    
    // artificial distance
    scalar_t d = dist_m;;
    
    map[mglevel->serialize(niv)] -= 2.0 * h / d * value;
    map[mglevel->serialize(n2iv)] += value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc0_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    // ABC0 == Dirichlet BC
    AmrIntVect_t niv = iv;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
    }
    
    map[mglevel->serialize(niv)] -= value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc1_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    // ABC1 == Robin BC (with normal derivatives) (i.e. Dirichlet BC + Neumann BC)
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    scalar_t sign = 1.0;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        n2iv[dir] = 1;
        
        sign = -1.0;
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    // correct cell size
    scalar_t h = mglevel->cellSize(dir);
    
    // coordinate of inner cell
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr) + sign * dist_m * h;
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr) + sign * dist_m * h;
#if AMREX_SPACEDIM == 3
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr) + sign * dist_m * h;
#endif
    scalar_t r = std::sqrt( AMREX_D_TERM(x * x, + y * y, + z * z) );
    
//     // artificial distance
//     scalar_t d = dist_m * nr[dir] * h; // + r;
    scalar_t d = r;
    
    map[mglevel->serialize(niv)] -= 2.0 * h / d * value;
    map[mglevel->serialize(n2iv)] += value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc2_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    AmrIntVect_t niv = iv;
    AmrIntVect_t n2iv = iv;
    
    scalar_t sign = 1.0;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        n2iv[dir] = 1;
        
        sign = -1.0;
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
        n2iv[dir] = nr[dir] - 2;
    }
    
    // correct cell size
    scalar_t h = mglevel->cellSize(dir);
    
    // coordinate of inner cell
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr) + sign * dist_m * h;
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr) + sign * dist_m * h;
#if AMREX_SPACEDIM == 3
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr) + sign * dist_m * h;
#endif
    scalar_t r = std::sqrt( AMREX_D_TERM(x * x, + y * y, + z * z) );
    
//     scalar_t d = dist_m * nr[dir] * h; // + r;
    scalar_t d = r;
    
    map[mglevel->serialize(niv)] += 2.0 * (d * d - h * h) / ( d * d + 2.0 * h * d ) * value;
    map[mglevel->serialize(n2iv)] += (2.0 * h - d) / (2.0 * h + d) * value;
}


template <class Level>
void AmrOpenBoundary<Level>::abc3_m(const AmrIntVect_t& iv,
                                    const lo_t& dir,
                                    umap_t& map,
                                    const scalar_t& value,
                                    Level* mglevel,
                                    const go_t* nr)
{
    AmrIntVect_t niv = iv;
    AmrIntVect_t n1iv = iv;
    AmrIntVect_t n2iv = iv;
    AmrIntVect_t n3iv = iv;
    AmrIntVect_t n4iv = iv;
    
    scalar_t sign = 1.0;
    
    if ( iv[dir] == -1 ) {
        // lower boundary
        niv[dir] = 0;
        n1iv[dir] = 1;
        n2iv[dir] = 2;
        n3iv[dir] = 3;
        n4iv[dir] = 4;
        
        
        
    } else {
        // upper boundary        
        niv[dir] = nr[dir] - 1;
        n1iv[dir] = nr[dir] - 2;
        n2iv[dir] = nr[dir] - 3;
        n3iv[dir] = nr[dir] - 4;
        n4iv[dir] = nr[dir] - 5;
        
        sign *= -1.0;
    }
    
    scalar_t h = mglevel->cellSize(dir);
    
//     scalar_t d = dist_m * nr[dir] * h; // + r;
    
    // coordinate of inner cell
    scalar_t x = this->coordinate_m(niv, 0, mglevel, nr);
    scalar_t y = this->coordinate_m(niv, 1, mglevel, nr);
#if AMREX_SPACEDIM == 3
    scalar_t z = this->coordinate_m(niv, 2, mglevel, nr);
#endif
    
    scalar_t d = std::sqrt( AMREX_D_TERM(x * x, + y * y, + z * z) );
    
    scalar_t hd = h + d;
    
    map[mglevel->serialize(niv)]  += (2.0 * d  / hd - 2.0 * h * h / (3.0 * d * hd) + sign * 5.0 * d * d / (18.0 * h * hd)) * value;
    map[mglevel->serialize(n1iv)] += (h / hd - d / hd - sign * d * d / ( hd * h) ) * value;
    map[mglevel->serialize(n2iv)] += sign * 4.0 * d * d / (3.0 * h * hd) * value;
    map[mglevel->serialize(n3iv)] -= sign * 7.0 * d * d / (9.0 * h * hd) * value;
    map[mglevel->serialize(n4iv)] += sign * d * d / (6.0 * h * hd) * value;
}


template <class Level>
typename AmrOpenBoundary<Level>::scalar_t
AmrOpenBoundary<Level>::coordinate_m(const AmrIntVect_t& iv,
                                     const lo_t& dir,
                                     Level* mglevel,
                                     const go_t* nr)
{
    /*
     * subtract / add half cell length in order to get
     * cell centered position
     */
    scalar_t h = mglevel->cellSize(dir);
    
    scalar_t lower = mglevel->geom.ProbLo(dir) + 0.5 * h;
    scalar_t upper = mglevel->geom.ProbHi(dir) - 0.5 * h;
    
    scalar_t m = (upper - lower) / ( nr[dir] - 1 );
    
    return m * iv[dir] + lower;
}

#endif
