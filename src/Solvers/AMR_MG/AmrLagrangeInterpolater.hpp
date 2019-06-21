#include "Utilities/OpalException.h"

#if AMREX_SPACEDIM == 3
template <class Level>
constexpr typename AmrLagrangeInterpolater<Level>::qpattern_t
    AmrLagrangeInterpolater<Level>::qpattern_ms;


template <class Level>
constexpr typename AmrLagrangeInterpolater<Level>::lpattern_t
    AmrLagrangeInterpolater<Level>::lpattern_ms;
#endif

//                                                      y_b   y_t
template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup1a_ms[2] = {0.25, -0.25};

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup2a_ms[2] = {0.75, 1.25};

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup1b_ms[2] = {-0.25, 0.25};

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup2b_ms[2] = {1.25, 0.75};

#if AMREX_SPACEDIM == 3
template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup3_ms[2] = {5.0 / 32.0, -3.0 / 32.0 };

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup3r_ms[2] = {-3.0 / 32.0, 5.0 / 32.0 };

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup4_ms[2] = {7.0 / 16.0, -9.0 / 16.0 };

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup4r_ms[2] = {-9.0 / 16.0, 7.0 / 16.0 };

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup5_ms[2] = {45.0 / 32.0, 21.0 / 32.0};

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup5r_ms[2] = {21.0 / 32.0, 45.0 / 32.0};

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::lookup6_ms = 15.0 / 16.0;

template <class Level>
const typename AmrLagrangeInterpolater<Level>::scalar_t
    AmrLagrangeInterpolater<Level>::factor_ms = 8.0 / 15.0;
#endif




template <class Level>
AmrLagrangeInterpolater<Level>::AmrLagrangeInterpolater(Order order)
    : AmrInterpolater<Level>( lo_t(order) + 1 )
{ }


template <class Level>
void AmrLagrangeInterpolater<Level>::stencil(
    const AmrIntVect_t& iv,
    const basefab_t& fab,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    Level* mglevel)
{
    
}


template <class Level>
void AmrLagrangeInterpolater<Level>::coarse(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    Level* mglevel)
{
    // polynomial degree = #points - 1
    switch ( this->nPoints_m - 1 ) {
        
        case Order::QUADRATIC:
            this->crseQuadratic_m(iv, map, scale, dir, shift, rfab, riv, mglevel);
            break;
        case Order::LINEAR:
            this->crseLinear_m(iv, map, scale, dir, shift, rfab, riv, mglevel);
            break;
        default:
            throw OpalException("AmrLagrangeInterpolater::coarse()",
                                "Not implemented interpolation");
    }
}


template <class Level>
void AmrLagrangeInterpolater<Level>::fine(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift,
    Level* mglevel)
{
    // polynomial degree = #points - 1
    switch ( this->nPoints_m - 1 ) {
        
        case Order::QUADRATIC:
            this->fineQuadratic_m(iv, map, scale, dir, shift, mglevel);
            break;
        case Order::LINEAR:
            this->fineLinear_m(iv, map, scale, dir, shift, mglevel);
            break;
        default:
            throw OpalException("AmrLagrangeInterpolater::fine()",
                                "Not implemented interpolation");
    }
}


template <class Level>
void AmrLagrangeInterpolater<Level>::fineLinear_m(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift,
    Level* mglevel)
{
    /*
     * computes the ghost cell directly
     */
    AmrIntVect_t tmp = iv;
    // first fine cell on refined coarse cell (closer to interface)
    tmp[dir] += shift;
    map[mglevel->serialize(tmp)] += 2.0 * scale;
    
    // second fine cell on refined coarse cell (further away from interface)
    tmp[dir] += shift;
    map[mglevel->serialize(tmp)] -= scale;
}


template <class Level>
void AmrLagrangeInterpolater<Level>::fineQuadratic_m(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift,
    Level* mglevel)
{
    AmrIntVect_t tmp = iv;
    // first fine cell on refined coarse cell (closer to interface)
    tmp[dir] += shift;
    map[mglevel->serialize(tmp)] += 2.0 / 3.0 * scale;
                        
    // second fine cell on refined coarse cell (further away from interface)
    tmp[dir] += shift;
    map[mglevel->serialize(tmp)] -= 0.2 * scale;
}


template <class Level>
void AmrLagrangeInterpolater<Level>::crseLinear_m(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    Level* mglevel)
{
#if AMREX_SPACEDIM == 2
    bool top = (riv[(dir+1)%AMREX_SPACEDIM] % 2 == 1);
    
    // right / upper / back
    AmrIntVect_t niv = iv;
    niv[(dir+1)%AMREX_SPACEDIM ] += 1;
    
    // left / lower / front
    AmrIntVect_t miv = iv;
    miv[(dir+1)%AMREX_SPACEDIM ] -= 1;
    
    // factor for fine
    scalar_t fac = 8.0 / 15.0 * scale;
    
    if ( mglevel->isValid(niv) && rfab(niv) != Level::Refined::YES ) {
        // check r / u / b --> 1: valid; 0: not valid
        map[mglevel->serialize(iv)] += fac * lookup2a_ms[top];
        map[mglevel->serialize(niv)] += fac * lookup1a_ms[top];
        
    } else if ( mglevel->isValid(miv) && rfab(miv) != Level::Refined::YES ) {
        // check l / f --> 1: valid; 0: not valid
        map[mglevel->serialize(iv)] += fac * lookup2a_ms[top];
        map[mglevel->serialize(miv)] += fac * lookup1a_ms[top];
        
    } else
        throw OpalException("AmrLagrangeInterpolater::crseLinear_m()",
                            "No valid interpolation scenario found!");
    
#elif AMREX_SPACEDIM == 3
    
    /*          x   y   z
     * ------------------
     * dir:     0   1   2
     * top1:    1   2   0   (i, L)
     * top2:    2   0   1   (j, K)
     */
    
    /* There are 4 coefficients from Lagrange interpolation.
     * Those are given by the product of one of
     * L0, L1 and one of K0, K1.
     * 
     * g(x, y) = f(x0, y0) * L0(x) * K0(y) +
     *           f(x0, y1) * L0(x) * K1(y) +
     *           f(x1, y0) * L1(x) * K0(y) +
     *           f(x1, y1) * L1(x) * K1(y) +
     */
    
    /*
     * check in 3x3 area (using iv as center) if 4 cells are not covered
     */
    lo_t d1 = (dir+1)%AMREX_SPACEDIM;
    lo_t d2 = (dir+2)%AMREX_SPACEDIM;
    
    lbits_t area;
    lo_t bit = 0;
    
    AmrIntVect_t tmp = iv;
    for (int i = -1; i < 2; ++i) {
        tmp[d1] += i;
        for (int j = -1; j < 2; ++j) {
            
            tmp[d2] += j;
            // make use of "short-circuit evaluation"
            area[bit] = ( mglevel->isValid(tmp) && rfab(tmp) );
            ++bit;
            
            // undo
            tmp[d2] -= j;
        }
        // undo
        tmp[d1] -= i;
    }
    
    qpattern_t::const_iterator pit = std::begin(this->lpattern_ms);
    
    while ( pit != std::end(this->lpattern_ms) ) {
        if ( *pit == (area & lbits_t(*pit)).to_ulong() )
            break;
        ++pit;
    }
    
    // factor for fine
    scalar_t fac = factor_ms * scale;
    
    scalar_t L[2] = {0.0, 0.0};
    lo_t top1 = riv[d1] % 2;
    
    scalar_t K[2] = {0.0, 0.0};
    lo_t top2 = riv[d2] % 2;
    
    lo_t begin[2] = { 0, 0 };
    lo_t end[2]   = { 0, 0 };
    
    switch ( *pit ) {
        case this->lpattern_ms[0]:
        {
            // corner top right pattern
            L[0] = lookup1a_ms[top1]; // L_{-1}
            L[1] = lookup2a_ms[top1]; // L_{0}
            begin[0] = -1;
            end[0]   =  0;
            
            K[0] = lookup1a_ms[top2]; // K_{-1}
            K[1] = lookup2a_ms[top2]; // K_{0}
            begin[1] = -1;
            end[1]   =  0;
            break;
        }
        case this->lpattern_ms[1]:
        {
            // corner bottom right pattern
            L[0] = lookup2b_ms[top1]; // L_{0}
            L[1] = lookup1b_ms[top1]; // L_{1}
            begin[0] = 0;
            end[0]   = 1;
            
            K[0] = lookup1a_ms[top2]; // K_{-1}
            K[1] = lookup2a_ms[top2]; // K_{0}
            begin[1] = -1;
            end[1]   =  0;
            break;
        }
        case this->lpattern_ms[2]:
        {
            // corner bottom left pattern
            L[0] = lookup2b_ms[top1]; // L_{0}
            L[1] = lookup1b_ms[top1]; // L_{1}
            begin[0] = 0;
            end[0]   = 1;
            
            K[0] = lookup2b_ms[top2]; // K_{0}
            K[1] = lookup1b_ms[top2]; // K_{1}
            begin[1] = 0;
            end[1]   = 1;
            break;
        }
        case this->lpattern_ms[3]:
        {
            // corner top left pattern
            L[0] = lookup1a_ms[top1]; // L_{-1}
            L[1] = lookup2a_ms[top1]; // L_{0}
            begin[0] = -1;
            end[0]   =  0;
            
            K[0] = lookup2b_ms[top2]; // K_{0}
            K[1] = lookup1b_ms[top2]; // K_{1}
            begin[1] = 0;
            end[1]   = 1;
            break;
        }
        default:
            throw OpalException("AmrLagrangeInterpolater::crseLinear_m()",
                                "No valid interpolation scenario found!");
    }
    
    /*
     * if pattern is known --> add stencil
     */
    AmrIntVect_t niv = iv;
    for (int i = begin[0]; i <= end[0]; ++i) {
        niv[d1] += i;
        for (int j = begin[1]; j <= end[1]; ++j) {
            niv[d2] += j;
            
            scalar_t value = fac * L[i-begin[0]] * K[j-begin[1]];
            if ( !mglevel->applyBoundary(niv, rfab, map, value) )
                map[mglevel->serialize(niv)] += value;
            
            // undo
            niv[d2] -= j;
        }
        // undo
        niv[d1] -= i;
    }
#else
    #error Lagrange interpolation: Only 2D and 3D are supported!
#endif
    // the neighbour cancels out
}


template <class Level>
void AmrLagrangeInterpolater<Level>::crseQuadratic_m(
    const AmrIntVect_t& iv,
    typename Level::umap_t& map,
    const typename Level::scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    Level* mglevel)
{
#if AMREX_SPACEDIM == 2
    
    bool top = (riv[(dir+1)%AMREX_SPACEDIM] % 2 == 1);
    
    // right / upper / back
    AmrIntVect_t niv = iv;
    niv[(dir+1)%AMREX_SPACEDIM ] += 1;
    
    // left / lower / front
    AmrIntVect_t miv = iv;
    miv[(dir+1)%AMREX_SPACEDIM ] -= 1;
    
    // 2nd right / upper / back
    AmrIntVect_t n2iv = niv;
    n2iv[(dir+1)%AMREX_SPACEDIM ] += 1;
    
    // 2nd left / lower / front
    AmrIntVect_t m2iv = miv;
    m2iv[(dir+1)%AMREX_SPACEDIM ] -= 1;
        
    /* 3 cases:
     * --------
     * r: right neighbour of iv (center of Laplacian)
     * u: upper neighbour of iv
     * b: back neighbour of iv
     * 
     * l: lower / left neighbour of iv
     * f: front neighbour of iv
     * 
     * -1 --> not valid, error happend
     * 0 --> r / u / b and 2nd r / u / b are valid
     * 1 --> direct neighbours (r / u / b and l / f) of iv are valid (has priority)
     * 2 --> l / f and 2nd l / f are valid
     */
    
    // check r / u / b --> 1: valid; 0: not valid
    bool rub = rfab(niv);
    
    // check l / f --> 1: valid; 0: not valid
    bool lf = rfab(miv);
    
    // check 2nd r / u / b
    bool rub2 = rfab(n2iv);
    
    // check 2nd l / f
    bool lf2 = rfab(m2iv);
    
    if ( rub && lf )
    {
        /*
         * standard case -1, +1 are not-refined nor at physical/mesh boundary
         */            
        // cell is not refined and not at physical boundary
        
        // y_t or y_b
        map[mglevel->serialize(iv)] += 0.5 * scale;
        
        //                             y_t          y_b
        scalar_t value = scale * ((top) ? 1.0 / 12.0 : -0.05);
        if ( !mglevel->applyBoundary(niv, rfab, map, value) )
            map[mglevel->serialize(niv)] += value;
        
        //                      y_t     y_b
        value = scale * ((top) ? -0.05 : 1.0 / 12.0);
        if ( !mglevel->applyBoundary(miv, rfab, map, value) )
            map[mglevel->serialize(miv)] += value;
        
    } else if ( rub && rub2 ) {
        /*
         * corner case --> right / upper / back + 2nd right / upper / back
         */
        //                     y_t          y_b
        scalar_t value = scale * ((top) ? 7.0 / 20.0 : 0.75);
        map[mglevel->serialize(iv)] += value;
        
        //                      y_t          y_b
        value = scale * ((top) ? 7.0 / 30.0 : -0.3);
        if ( !mglevel->applyBoundary(niv, rfab, map, value) )
            map[mglevel->serialize(niv)] += value;
        
        //                      y_t     y_b
        value = scale * ((top) ? -0.05 : 1.0 / 12.0);
        if ( !mglevel->applyBoundary(n2iv, rfab, map, value) )
            map[mglevel->serialize(n2iv)] += value;
        
    } else if ( lf && lf2 ) {
        /*
         * corner case --> left / lower / front + 2nd left / lower / front
         */
        //                             y_t    y_b
        scalar_t value = scale * ((top) ? 0.75 : 7.0 / 20.0);
        map[mglevel->serialize(iv)] += value;
        
        //                      y_t           y_b
        value = scale * ((top) ? -0.3 :  7.0 / 30);
        if ( !mglevel->applyBoundary(miv, rfab, map, value) )
            map[mglevel->serialize(miv)] += value;
        
        //                      y_t          y_b
        value = scale * ((top) ? 1.0 / 12.0 : -0.05);
        if ( !mglevel->applyBoundary(m2iv, rfab, map, value) )
            map[mglevel->serialize(m2iv)] += value;
        
    } else {
        /* last trial: linear Lagrange interpolation
         * --> it throws an error if not possible
         */
        this->crseLinear_m(iv, map, scale, dir, shift, rfab, riv, mglevel);
    }
    
#elif AMREX_SPACEDIM == 3
    
    /*          x   y   z
     * ------------------
     * dir:     0   1   2
     * top1:    1   2   0   (i, L)
     * top2:    2   0   1   (j, K)
     */
    
    /* There are 9 coefficients from Lagrange interpolation.
     * Those are given by the product of one of
     * L0, L1, L2 and one of K0, K1, K2.
     * 
     * g(x, y) = f(x0, y0) * L0(x) * K0(y) +
     *           f(x0, y1) * L0(x) * K1(y) +
     *           f(x0, y2) * L0(x) * K2(y) +
     *           f(x1, y0) * L1(x) * K0(y) +
     *           f(x1, y1) * L1(x) * K1(y) +
     *           f(x1, y2) * L1(x) * K2(y) +
     *           f(x2, y0) * L2(x) * K0(y) +
     *           f(x2, y1) * L2(x) * K1(y) +
     *           f(x2, y2) * L2(x) * K2(y) +
     */
    
    
    /*
     * check in 5x5 area (using iv as center) if 9 cells are not covered
     */
    lo_t d1 = (dir+1)%AMREX_SPACEDIM;
    lo_t d2 = (dir+2)%AMREX_SPACEDIM;
    
    qbits_t area;
    lo_t bit = 0;
    
    AmrIntVect_t tmp = iv;
    for (int i = -2; i < 3; ++i) {
        tmp[d1] += i;
        for (int j = -2; j < 3; ++j) {
            
            tmp[d2] += j;
            
            area[bit] = ( mglevel->isValid(tmp) && rfab(tmp) );
            ++bit;
            
            // undo
            tmp[d2] -= j;
        }
        // undo
        tmp[d1] -= i;
    }
    
    qpattern_t::const_iterator pit = std::begin(this->qpattern_ms);

    while ( pit != std::end(this->qpattern_ms) ) {
        if ( *pit == (area & qbits_t(*pit)).to_ulong() )
            break;
        ++pit;
    }
    
    // factor for fine
    scalar_t fac = factor_ms * scale;
    
    scalar_t L[3] = {0.0, 0.0, 0.0};
    lo_t top1 = riv[d1] % 2;
    
    scalar_t K[3] = {0.0, 0.0, 0.0};
    lo_t top2 = riv[d2] % 2;
    
    lo_t begin[2] = { 0, 0 };
    lo_t end[2]   = { 0, 0 };
    
    switch ( *pit ) {
        case this->qpattern_ms[0]:
        {
            // cross pattern
            L[0] = lookup3_ms[top1];  // L_{-1}
            L[1] = lookup6_ms;        // L_{0}
            L[2] = lookup3r_ms[top1]; // L_{1}
            begin[0] = -1;
            end[0]   =  1;
            
            K[0] = lookup3_ms[top2];  // K_{-1}
            K[1] = lookup6_ms;        // K_{0}
            K[2] = lookup3r_ms[top2]; // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            break;
        }
        case this->qpattern_ms[1]:
        {
            // T pattern
            L[0] = lookup3r_ms[top1]; // L_{-2}
            L[1] = lookup4_ms[top1];  // L_{-1}
            L[2] = lookup5r_ms[top1]; // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = lookup3_ms[top2];  // K_{-1}
            K[1] = lookup6_ms;        // K_{0}
            K[2] = lookup3r_ms[top2]; // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            break;
        }
        case this->qpattern_ms[2]:
        {
            // right hammer pattern
            L[0] = lookup3_ms[top1];  // L_{-1}
            L[1] = lookup6_ms;        // L_{0}
            L[2] = lookup3r_ms[top1]; // L_{1}
            begin[0] = -1;
            end[0] = 1;
            
            K[0] = lookup3r_ms[top2]; // K_{-2}
            K[1] = lookup4_ms[top2];  // K_{-1}
            K[2] = lookup5r_ms[top2]; // K_{0}
            begin[1] = -2;
            end[1]   = 0;
            break;
        }
        case this->qpattern_ms[3]:
        {
            // T on head pattern
            L[0] = lookup5_ms[top1];  // L_{0}
            L[1] = lookup4r_ms[top1]; // L_{1}
            L[2] = lookup3_ms[top1];  // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = lookup3_ms[top2];  // K_{-1}
            K[1] = lookup6_ms;        // K_{0}
            K[2] = lookup3r_ms[top2]; // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            break;
        }
        case this->qpattern_ms[4]:
        {
            // left hammer pattern
            L[0] = lookup3_ms[top1];  // L_{-1}
            L[1] = lookup6_ms;        // L_{0}
            L[2] = lookup3r_ms[top1]; // L_{1}
            begin[0] = -1;
            end[0] = 1;
            
            K[0] = lookup5_ms[top2];  // K_{0}
            K[1] = lookup4r_ms[top2]; // K_{1}
            K[2] = lookup3_ms[top2];  // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[5]:
        {
            // upper left corner pattern
            L[0] = lookup3r_ms[top1]; // L_{-2}
            L[1] = lookup4_ms[top1];  // L_{-1}
            L[2] = lookup5r_ms[top1]; // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = lookup5_ms[top2];  // K_{0}
            K[1] = lookup4r_ms[top2]; // K_{1}
            K[2] = lookup3_ms[top2];  // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[6]:
        {
            // upper right corner pattern
            L[0] = lookup3r_ms[top1]; // L_{-2}
            L[1] = lookup4_ms[top1];  // L_{-1}
            L[2] = lookup5r_ms[top1]; // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = lookup3r_ms[top2]; // K_{-2}
            K[1] = lookup4_ms[top2];  // K_{-1}
            K[2] = lookup5r_ms[top2]; // K_{0}
            begin[1] = -2;
            end[1]   =  0;
            break;
        }
        case this->qpattern_ms[7]:
        {
            // mirrored L pattern
            L[0] = lookup5_ms[top1];  // L_{0}
            L[1] = lookup4r_ms[top1]; // L_{1}
            L[2] = lookup3_ms[top1];  // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = lookup3r_ms[top2]; // K_{-2}
            K[1] = lookup4_ms[top2];  // K_{-1}
            K[2] = lookup5r_ms[top2]; // K_{0}
            begin[1] = -2;
            end[1]   =  0;
            break;
        }
        case this->qpattern_ms[8]:
        {
            // L pattern
            L[0] = lookup5_ms[top1];  // L_{0}
            L[1] = lookup4r_ms[top1]; // L_{1}
            L[2] = lookup3_ms[top1];  // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = lookup5_ms[top2];  // K_{0}
            K[1] = lookup4r_ms[top2]; // K_{1}
            K[2] = lookup3_ms[top2];  // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        default:
        {
            /* unknown pattern --> last trial: linear Lagrange interpolation
             * --> it throws an error if not possible
             */
            this->crseLinear_m(iv, map, scale, dir, shift, rfab, riv, mglevel);
            return;
        }
    }
    
    /*
     * if pattern is known --> add stencil
     */
    AmrIntVect_t tmp1 = iv;
    for (int i = begin[0]; i <= end[0]; ++i) {
        tmp1[d1] += i;
        for (int j = begin[1]; j <= end[1]; ++j) {
            tmp1[d2] += j;
            
            scalar_t value = fac * L[i-begin[0]] * K[j-begin[1]];
            if ( !mglevel->applyBoundary(tmp1, rfab, map, value) )
                map[mglevel->serialize(tmp1)] += value;
            
            // undo
            tmp1[d2] -= j;
        }
        // undo
        tmp1[d1] -= i;
    }
    
#else
    #error Lagrange interpolation: Only 2D and 3D are supported!
#endif
}
