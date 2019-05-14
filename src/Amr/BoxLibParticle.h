#ifndef BOXLIB_PARTICLE_H
#define BOXLIB_PARTICLE_H

#include "AmrParticle/AmrParticleBase.h"

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>

/*!
 * Particle class for AMReX. It works together with BoxLibLayout.
 * 
 * Particle class that does the scatter and gather operations of attributes
 * to and from grid. In contrast to Ippl where it is implemented in the
 * attribute class, i.e. src/ippl/src/Particle/ParticleAttrib.h we do it in
 * the particle class. This way Ippl does not need to be modified if another
 * AMR framework is used.
 */
template<class PLayout>
class BoxLibParticle : public virtual AmrParticleBase<PLayout>
{
public:
    typedef typename AmrParticleBase<PLayout>::ParticlePos_t             ParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::ParticleIndex_t           ParticleIndex_t;
    typedef typename AmrParticleBase<PLayout>::SingleParticlePos_t       SingleParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::AmrField_t                AmrField_t;
    // Array<std::unique_ptr<MultiFab> >
    typedef typename AmrParticleBase<PLayout>::AmrScalarFieldContainer_t AmrScalarFieldContainer_t;
    typedef typename AmrParticleBase<PLayout>::AmrVectorFieldContainer_t AmrVectorFieldContainer_t;
    typedef typename AmrParticleBase<PLayout>::AmrVectorField_t          AmrVectorField_t;
    typedef typename AmrParticleBase<PLayout>::ParticleLevelCounter_t    ParticleLevelCounter_t;
    
    typedef typename PLayout::AmrProcMap_t  AmrProcMap_t;
    typedef typename PLayout::AmrGrid_t     AmrGrid_t;
    typedef typename PLayout::AmrGeometry_t AmrGeometry_t;
    typedef typename PLayout::AmrIntVect_t  AmrIntVect_t;
    typedef typename PLayout::AmrBox_t      AmrBox_t;
    typedef typename PLayout::AmrReal_t     AmrReal_t;
    
    typedef amrex::FArrayBox                FArrayBox_t;
    
public:
    BoxLibParticle();
    
    /*!
     * @param layout that does the particle-to-core management
     */
    BoxLibParticle(PLayout *layout);
    
    /*!
     * Multi-level scatter.
     * Scatter the data from the given attribute onto the given field, using
     * the given position attribute. It calls the AMReX method.
     * 
     * @param attrib to scatter onto grid
     * @param f field on grid
     * @param pp particle position (not used for AMReX call)
     * @param lbase base level we want to start
     * @param lfine finest level we want to stop
     * @param pbin the particle bin attribute
     * @param bin to scatter (default: -1 --> scatter all particles)
     */
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrScalarFieldContainer_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int lbase, int lfine,
                 const ParticleAttrib<int>& pbin, int bin = -1);
    
    
    /*!
     * Single-level scatter.
     * Scatter the data from the given attribute onto the given field, using
     * the given position attribute. It calls the AMReX methods.
     * 
     * @param attrib to scatter onto grid
     * @param f field on grid
     * @param pp particle position (not used for AMReX call)
     * @param pbin the particle bin attribute
     * @param bin to scatter (default: -1 --> scatter all particles)
     * @param level for which we put particles onto the grid
     */
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrField_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 const ParticleAttrib<int>& pbin, int bin = -1,
                 int level = 0);
    
    /*!
     * Multi-level gather.
     * Gather the data from the given Field into the given attribute, using
     * the given Position attribute.
     * 
     * @param attrib to gather from grid
     * @param f vector field on grid
     * @param pp particle position (not used for AMReX call)
     * @param lbase base level to gather from
     * @param lfine finest level to gather from
     */
    template <class FT, unsigned Dim, class PT>
    void gather(ParticleAttrib<FT>& attrib, AmrVectorFieldContainer_t& f,
                ParticleAttrib<Vektor<PT, Dim> >& pp,
                int lbase, int lfine);
    
private:
    /*
     * AMReX functions adjusted to work with Ippl
     */
    
    /*!
     * Multi-level scatter (adjusted from AMReX).
     * 
     * @param pa is the attribute to scatter onto the grid
     * @param mf_to_be_filled is the MultiFab container to be filled
     * (i.e. grid data)
     * @param lev_min level we want to start
     * @param ncomp is the number of components of MultiFab (equal to 1)
     * @param finest_level level we want to end
     */
    template <class AType>
    void AssignDensityFort(ParticleAttrib<AType> &pa,
                           AmrScalarFieldContainer_t& mf_to_be_filled, 
                           int lev_min, int ncomp, int finest_level,
                           const ParticleAttrib<int>& pbin, int bin = -1) const;
    
    /*!
     * Multi-level gather (adjusted from AMReX).
     * 
     * @param pa is the attribute to gather to.
     * @param mesh_data where the information is
     * @param lev_min level to start
     * @param lev_max level to end
     */
    template <class AType>
    void InterpolateFort(ParticleAttrib<AType> &pa,
                         AmrVectorFieldContainer_t& mesh_data, 
                         int lev_min, int lev_max);
    
    /*!
     * Single-level gather (adjusted from AMReX).
     * 
     * @param pa is the attribute to be updated
     * @param mesh_data where the information is taken from
     * @param lev for which we get the mesh data
     */
    template <class AType>
    void InterpolateSingleLevelFort(ParticleAttrib<AType> &pa, AmrVectorField_t& mesh_data, int lev);
    
    
    /*!
     * Multi-level gather.
     * 
     * @param pa is the attribute to be updated
     * @param mesh_data where the information is taken from
     * @param lev for which we get the mesh data
     */
    template <class AType>
    void InterpolateMultiLevelFort(ParticleAttrib<AType> &pa,
                                   AmrVectorFieldContainer_t& mesh_data,
                                   int lev);
    
    
    /*!
     * Single-level scatter (adjusted from AMReX).
     * 
     * @param pa is the attribute to scatter onto the grid
     * @param mf where attribute is scatterd to
     * @param level where we want to scatter
     * @param ncomp is the number of the component in the MultiFab (ncomp = 1)
     * @param particle_lvl_offset is zero
     */
    template <class AType>
    void AssignCellDensitySingleLevelFort(ParticleAttrib<AType> &pa, AmrField_t& mf, int level,
                                          const ParticleAttrib<int>& pbin, int bin = -1,
                                          int ncomp=1, int particle_lvl_offset = 0) const;
    
private:
    IpplTimings::TimerRef AssignDensityTimer_m;
};


#include "BoxLibParticle.hpp"

#endif
