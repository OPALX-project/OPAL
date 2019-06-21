#ifndef AMR_LAGRANGE_INTERPOLATER_H
#define AMR_LAGRANGE_INTERPOLATER_H

#include "AmrInterpolater.h"

#if AMREX_SPACEDIM == 3
    #include <bitset>
    #include <iterator>
    #include <array>
#endif

#include "Ippl.h"

template <class Level>
class AmrLagrangeInterpolater : public AmrInterpolater<Level>
{
public:
    
    typedef typename Level::go_t        go_t;
    typedef typename Level::lo_t        lo_t;
    typedef typename Level::scalar_t    scalar_t;
    typedef typename Level::umap_t      umap_t;
    typedef typename Level::basefab_t   basefab_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;

    enum Order {
        LINEAR = 1,
        QUADRATIC
    };
    
#if AMREX_SPACEDIM == 3
    typedef std::bitset<25> qbits_t; ///< for checking the neighbour cells (quadratic)
    typedef std::bitset<9> lbits_t; ///< for checking the neighbour cells (linear)
    typedef std::array<unsigned int long, 9> qpattern_t;    ///< quadratic pattern
    typedef std::array<unsigned int long, 4> lpattern_t;    ///< linear pattern
#endif
    
public:
    
    AmrLagrangeInterpolater(Order order);
    
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
              lo_t dir, lo_t shift,
              Level* mglevel);
    
private:
    
    /*!
     * First order interpolation on fine cell interface side
     * 
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
    void fineLinear_m(const AmrIntVect_t& iv,
                      umap_t& map,
                      const scalar_t& scale,
                      lo_t dir, lo_t shift,
                      Level* mglevel);
    
    /*!
     * Second order interpolation on fine cell interface side
     * 
     * @param iv is the fine ghost cell at the interface (on coarse cell that is not
     * refined)
     * @param map with global matrix indices of fine level cells and
     * values matrix entries of fine level cells (coefficients)
     * @param scale of matrix values
     * @param dir direction of interface (0 "horizontal", 1 "vertical", 2 "longitudinal")
     * @param shift is either -1 or 1. If the refined coarse cell is on the left / lower / front
     * side, shift is equal to -1, otherwise the interface is on the right / upper / back side
     * and the value is 1.
     * @param mglevel used to get the global indices and refinement ratio among levels,
     * and boundary avlues at physical domain, e.g. Dirichlet, open BC
     */
    void fineQuadratic_m(const AmrIntVect_t& iv,
                         umap_t& map,
                         const scalar_t& scale,
                         lo_t dir, lo_t shift,
                         Level* mglevel);
    
    /*!
     * First oder interpolation on coarse cell interface side
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * values matrix entries of coarse level cells (coefficients)
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
    void crseLinear_m(const AmrIntVect_t& iv,
                      umap_t& map,
                      const scalar_t& scale,
                      lo_t dir, lo_t shift, const basefab_t& rfab,
                      const AmrIntVect_t& riv,
                      Level* mglevel);
    
    /*!
     * Second order interpolation on coarse cell interface side
     * 
     * @param iv is the coarse cell at the interface (center cell of Laplacian)
     * @param map with global matrix indices of coarse level cells and
     * values matrix entries of coarse level cells (coefficients)
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
    void crseQuadratic_m(const AmrIntVect_t& iv,
                         umap_t& map,
                         const scalar_t& scale,
                         lo_t dir, lo_t shift, const basefab_t& rfab,
                         const AmrIntVect_t& riv,
                         Level* mglevel);
    
private:
#if AMREX_SPACEDIM == 3
    static constexpr qpattern_t qpattern_ms {
        473536,                             ///< cross pattern (case 0)
        14798,                              ///< T pattern (case 1)
        236768,                             ///< right hammer pattern (case 2)
        15153152,                           ///< T on head pattern (case 3)
        947072,                             ///< left hammer pattern (case 4)
        29596,                              ///< upper left corner pattern (case 5)
        7399,                               ///< upper right corner pattern (case 6)
        7576576,                            ///< mirrored L pattern (case 7)
        30306304                            ///< L pattern (case 8)
    };
    
    static constexpr lpattern_t lpattern_ms {
        27,                                 ///< corner top right pattern
        216,                                ///< corner bottom right pattern
        432,                                ///< corner bottom left pattern
        54                                  ///< corner top left pattern
    };
#endif

    // y_b   y_t
    static const scalar_t lookup1a_ms[2];
    static const scalar_t lookup2a_ms[2];
    static const scalar_t lookup1b_ms[2];
    static const scalar_t lookup2b_ms[2];
#if AMREX_SPACEDIM == 3
    static const scalar_t lookup3_ms[2];
    static const scalar_t lookup3r_ms[2];
    static const scalar_t lookup4_ms[2];
    static const scalar_t lookup4r_ms[2];
    static const scalar_t lookup5_ms[2];
    static const scalar_t lookup5r_ms[2];
    static const scalar_t lookup6_ms;
    static const scalar_t factor_ms;
#endif
};

#include "AmrLagrangeInterpolater.hpp"

#endif
