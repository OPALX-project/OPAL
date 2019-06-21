#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"
#include "AmrParticle/AmrParticleBase.h"

#include "Amr/AmrDefs.h"

#include <AMReX_ParGDB.H>

/*!
 * In contrast to AMReX, OPAL is optimized for
 * distribution particles to cores. In AMReX the ParGDB object
 * is responsible for the particle to core distribution. This
 * layout is derived from this object and does all important
 * bunch updates. It is the interface for AMReX and Ippl.
 * 
 * In AMReX the geometry, i.e. physical domain, is fixed
 * during the whole computation. Particles leaving the domain
 * would be deleted. In order to prevent this we map the particles
 * onto the domain \f$[-1, 1]^3\f$. Furthermore, it makes sure
 * that we have enougth grid points to represent the bunch
 * when its charges are scattered on the grid for the self-field
 * computation.
 * 
 * The self-field computation and the particle-to-core update
 * are performed in the particle mapped domain.
 */
template<class T, unsigned Dim>
class BoxLibLayout : public ParticleAmrLayout<T, Dim>,
                     public amrex::ParGDB
{
    
public:
    typedef typename ParticleAmrLayout<T, Dim>::pair_t pair_t;
    typedef typename ParticleAmrLayout<T, Dim>::pair_iterator pair_iterator;
    typedef typename ParticleAmrLayout<T, Dim>::SingleParticlePos_t SingleParticlePos_t;
    typedef typename ParticleAmrLayout<T, Dim>::Index_t Index_t;
    
    typedef amr::AmrField_t                AmrField_t;
    typedef amr::AmrVectorField_t          AmrVectorField_t;
    typedef amr::AmrScalarFieldContainer_t AmrScalarFieldContainer_t;
    typedef amr::AmrVectorFieldContainer_t AmrVectorFieldContainer_t;
    typedef typename ParticleAmrLayout<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef ParticleAttrib<Index_t> ParticleIndex_t;
    
    typedef amr::AmrProcMap_t           AmrProcMap_t;
    typedef amr::AmrGrid_t              AmrGrid_t;
    typedef amr::AmrGeometry_t          AmrGeometry_t;
    typedef amr::AmrGeomContainer_t     AmrGeomContainer_t;
    typedef amr::AmrGridContainer_t     AmrGridContainer_t;
    typedef amr::AmrProcMapContainer_t  AmrProcMapContainer_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    typedef amr::AmrIntVectContainer_t  AmrIntVectContainer_t;
    typedef amr::AmrIntArray_t          AmrIntArray_t;
    typedef amr::AmrDomain_t            AmrDomain_t;
    typedef amr::AmrBox_t               AmrBox_t;
    typedef amr::AmrReal_t              AmrReal_t;
    
    typedef amrex::BaseFab<int>         basefab_t;
    typedef amrex::FabArray<basefab_t>  mask_t;
    
    /*!
     * Lower physical domain boundary (each dimension). It has to be
     * smaller than -1 since all particles are within \f$[-1, 1]^3\f$.
     * The real computational domain is multiplied with the mesh
     * enlargement factor (in [%]) in BoxLibLayout::initBaseBox_m().
     */
    static Vector_t lowerBound;
    
    /*! Upper physical domain boundary (each dimension). It has to be
     * greater than 1 since all particles are within \f$[-1, 1]^3\f$.
     * The real computational domain is multiplied with the mesh
     * enlargement factor (in [%]) in BoxLibLayout::initBaseBox_m().
     */
    static Vector_t upperBound;

public:
    
    /*!
     * Initializes default Geometry, DistributionMapping and BoxArray.
     */
    BoxLibLayout();
    
    /*!
     * Given a layout it copies that.
     */
    BoxLibLayout(const BoxLibLayout* layout_p);

    /*!
     * @param nGridPoints per dimension (nx, ny, nz / nt)
     * @param maxGridSize for all levels.
     */
    BoxLibLayout(int nGridPoints, int maxGridSize);
    
    /*!
     * Single-level constructor.
     * 
     * @param geom specifies the box domain
     * @param dmap is the distribution map for grids
     * @param ba is the array of boxes for a level
     */
    BoxLibLayout(const AmrGeometry_t &geom,
                 const AmrProcMap_t &dmap,
                 const AmrGrid_t &ba);
    
    /*!
     * Multi-level constructor.
     * 
     * @param geom is basically the physical domain storing
     * the mesh spacing per level
     * @param dmap are all distribution maps of grids to core
     * @param ba are all boxes of all levels
     * @param rr is the refinement ratio among the levels
     * (always the ratio from l to l+1)
     */
    BoxLibLayout(const AmrGeomContainer_t &geom,
                 const AmrProcMapContainer_t &dmap,
                 const AmrGridContainer_t &ba,
                 const AmrIntArray_t &rr);
    
    
    /*
     * Overloaded functions of ParticleAmrLayout
     */
    
    /*!
     * This method is used when creating the AMR object. OPAL
     * takes the input argument BBOXINCR that is specified in
     * the field solver command. Up to this point the
     * AMR object is not yet initialized. After that this
     * method shouldn't be called anymore.
     * @param dh is the mesh enlargement factor
     */
    void setBoundingBox(double dh);
    
    /*!
     * The Poisson computation domain is per default [-1,1]^3.
     * With this method this can be changed in order to account
     * for the different extent per direction.
     * @param ratio, if e.g. ratio = [1, 2, 0.75], then the box is
     * [-1, 1] x [-2, 2] x [-0.75, 0.75]
     */
    void setDomainRatio(const std::vector<double>& ratio);


    /*
     * Functions of IpplParticleBase
     */
    
    /*!
     * This method shouldn't be called. Otherwise
     * it throws an exception.
     */
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData,
                const ParticleAttrib<char>* canSwap = 0);
    
    /*!
     * The proper update method for AMR.
     * 
     * @param PData is basically the bunch
     * @param lev_min base level to update
     * @param lev_max finest level to update
     * (if -1 update all levels starting from lev_min)
     * @param isRegrid in regrid mode a level might be deleted, thus
     * the particle level counter might be invalidated --> we need to
     * iterate through all particles.
     */
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData,
                int lev_min = 0, int lev_max = -1, bool isRegrid = false);
    
    
    /*
     * Functions from AMReX that are adjusted to work with Ippl AmrParticleBase class
     */
    
    /*!
     * Get the cell of a particle
     * 
     * @param p is the particle data
     * @param ip is the local index of the particle in the container
     * @param level of the particle
     */
    AmrIntVect_t Index (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                        const unsigned int ip, int level) const;
    
    /*!
     * Get the cell of a particle
     * 
     * @param R is the position of a particle
     * @param lev is the level
     */
    AmrIntVect_t Index (SingleParticlePos_t &R, int lev) const;
    
    
    /*
     * Additional methods
     */
    
    
    /*!
     * Build mask for a level used for interpolation from
     * grid to particles to reduce spurious self field
     * forces near coarse-fine interfaces.
     */
    void buildLevelMask(int lev, const int ncells = 1);

    void clearLevelMask(int lev);

    const std::unique_ptr<mask_t>& getLevelMask(int lev) const;
    
    
    /*!
     * The particles live initially on the coarsest level.
     * Furthermore, the order the OPAL input file is parsed
     * does not allow us to know the maximum level of the
     * computation. This is known after that the FieldSolver
     * is initialized. Therefore, we need to update the
     * size of the ParGDB containers.
     * 
     * @param maxLevel is set when the FieldSolver is
     * initialized
     */
    void resize(int maxLevel) {
        int length = maxLevel + 1;
        this->m_geom.resize(length);
        this->m_dmap.resize(length);
        this->m_ba.resize(length);
        this->m_nlevels = length;
        this->refRatio_m.resize(maxLevel);
//         this->m_rr.resize(maxLevel);
        this->maxLevel_m = maxLevel;
    }
    

    /*!
     * Set the geometry of the problem. It is called in
     * AmrBoxLib::initBaseLevel_m().
     * 
     * @param geom geometry of all levels
     */
    void define(const AmrGeomContainer_t& geom) {
        for (unsigned int i = 0; i < geom.size(); ++i)
            this->m_geom[i] = geom[i];
    }
    
    
    /*!
     * Set the refinement ratios. It is called in
     * AmrBoxLib::initBaseLevel_m().
     * 
     * @param refRatio among levels
     */
    void define(const AmrIntVectContainer_t& refRatio) {
        for (unsigned int i = 0; i < refRatio.size(); ++i) {
            refRatio_m[i] = refRatio[i];
        }
    }
    
    /*
     * ParGDB overwritten functions
     */
    
    /*!
     * Check if an AMR level is well defined
     * 
     * @param level to check
     */
    inline bool LevelDefined (int level) const;
    
    /*!
     * @returns the current finest level
     */
    inline int finestLevel () const;
    
    /*!
     * @returns the maximum level of simulation
     */
    inline int maxLevel () const;
    
    /*!
     * @param level
     * @returns the refinement ratio of this level to the next
     * higher one
     */
    inline AmrIntVect_t refRatio (int level) const;
    
    /*!
     * @param level
     * @returns the maximum refinement ratio among all directions
     * for the given level.
     */
    inline int MaxRefRatio (int level) const;
    
private:
    
    /*!
     * Set up the box for the whole computation.
     * The AMR object owning the bunch is not yet initialized.
     * 
     * @param nGridPoints per dimension (nx, ny, nz / nt)
     * @param maxGridSize for all levels
     * @param dh is the mesh enlargement factor
     */
    void initBaseBox_m(int nGridPoints, int maxGridSize, double dh = 0.04);
    
    
    /*
     * Functions from AMReX that are adjusted to work with Ippl AmrParticleBase class
     */
    
    /*!
     * Function from AMReX adjusted to work with Ippl AmrParticleBase class
     * Checks/sets a particles location on levels lev_min and higher.
     * 
     * @param p is the bunch information
     * @param ip is the local (i.e. to a core) particle index
     * @param lev_min to check
     * @param lev_max to check
     * @param nGrow is the number of ghost cells
     * @returns false if the particle does not exist on that level.
     */
    bool Where (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                const unsigned int ip, 
                int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

    /*!
     * Function from AMReX adjusted to work with Ippl AmrParticleBase class
     * Checks/sets whether the particle has crossed a periodic boundary in such a way
     * that it is on levels lev_min and higher.
     * 
     * @param prt is the bunch information
     * @param ip is the local (i.e. to a core) particle index
     * @param lev_min to check
     * @param lev_max to check
     * @returns true if mapped to the other side.
     */
    bool EnforcePeriodicWhere (AmrParticleBase< BoxLibLayout<T,Dim> >& prt,
                               const unsigned int ip,
                               int lev_min = 0, int lev_max = -1) const;

    /*!
     * Function from AMReX adjusted to work with Ippl AmrParticleBase class
     * 
     * Move the particle to the opposite side of the domain
     * 
     * @param R is the particle position
     * @returns true if the particle was shifted.
     */
    bool PeriodicShift (SingleParticlePos_t R) const;
    
    /*!
     * Function from AMReX adjusted to work with Ippl AmrParticleBase class
     * 
     * @param p is basically the bunch
     * @param ip is the local particle index
     * @param lev_min to check
     * @param lev_max to check
     * @param nGrow is the number of ghost cells
     */
    void locateParticle(AmrParticleBase< BoxLibLayout<T,Dim> >& p, 
                        const unsigned int ip,
                        int lev_min, int lev_max, int nGrow) const;
    
private:
    
    // don't use m_rr from ParGDB since it is the same refinement in all directions
    AmrIntVectContainer_t refRatio_m;   /// Refinement ratios [0:finest_level-1]
    
    /* mask to reduce spurious self-field forces at
     * coarse-fine interfaces
     */
    std::vector<std::unique_ptr<mask_t> > masks_m;
};

#include "BoxLibLayout.hpp"

#endif
