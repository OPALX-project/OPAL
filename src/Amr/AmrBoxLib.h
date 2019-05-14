#ifndef AMR_BOXLIB_H
#define AMR_BOXLIB_H

#include "Amr/AmrObject.h"

class AmrPartBunch;

// AMReX headers
#include <AMReX_AmrMesh.H>
#include <AMReX.H>

/*!
 * Concrete AMR object. It is based on the
 * <a href="https://ccse.lbl.gov/AMReX/">AMReX</a> library
 * developed at LBNL. This library is the successor of
 * <a href="https://ccse.lbl.gov/BoxLib/">BoxLib</a>.
 */
class AmrBoxLib : public AmrObject,
                  public amrex::AmrMesh
{
    
public:
    typedef amr::AmrField_t                 AmrField_t;
    typedef amr::AmrScalarFieldContainer_t  AmrScalarFieldContainer_t;
    typedef amr::AmrVectorFieldContainer_t  AmrVectorFieldContainer_t;
    typedef amr::AmrGeomContainer_t         AmrGeomContainer_t;
    typedef amr::AmrGridContainer_t         AmrGridContainer_t;
    typedef amr::AmrProcMapContainer_t      AmrProcMapContainer_t;
    typedef amr::AmrDomain_t                AmrDomain_t;
    typedef amr::AmrIntArray_t              AmrIntArray_t;
    typedef amr::AmrReal_t                  AmrReal_t;
    typedef amr::AmrGrid_t                  AmrGrid_t;
    typedef amr::AmrProcMap_t               AmrProcMap_t;
    typedef amr::AmrGeometry_t              AmrGeometry_t;
    typedef amr::AmrIntVect_t               AmrIntVect_t;
    
    typedef amrex::FArrayBox                FArrayBox_t;
    typedef amrex::Box                      Box_t;
    typedef amrex::TagBox                   TagBox_t;
    typedef amrex::TagBoxArray              TagBoxArray_t;
    typedef amrex::MFIter                   MFIter_t;
    
    
    /*!
     * Used for the redistribution of grids
     */
    enum Strategy {
        RANK_ZERO = 0,      // all grids to processor zero
        PFC       = 1,
        RANDOM    = 2,
        KNAPSACK  = 3
    };
    
public:
    
    /*!
     * See other constructors documentation for further info.
     * @param domain is the physical domain of the problem
     * @param nGridPts per dimension (nx, ny, nz / nt)
     * @param maxLevel of mesh refinement
     * @param bunch is used when we tag for charges per cell
     */
    AmrBoxLib(const AmrDomain_t& domain,
              const AmrIntArray_t& nGridPts,
              int maxLevel,
              AmrPartBunch* bunch_p);
    
    /*!
     * Create a new object
     * @param info are the initial informations to construct an AmrBoxLib object.
     * It is set in src/Solvers/FieldSolver.cpp
     * @param bunch_p pointing to. Is used for getting the geometry (i.e. physical
     * domain of the problem).
     */
    static std::unique_ptr<AmrBoxLib> create(const AmrInfo& info,
                                             AmrPartBunch* bunch_p);
    
    /*!
     * Inherited from AmrObject
     */
    void regrid(double time);

    void getGridStatistics(std::map<int, long>& gridPtsPerCore,
                           std::vector<int>& gridsPerLevel) const;
    
    /*!
     * Initial gridding. Sets up all levels.
     */
    void initFineLevels();

    VectorPair_t getEExtrema();
    
    double getRho(int x, int y, int z);
    
    void computeSelfFields();
    
    void computeSelfFields(int bin);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int bin);
    
    void updateMesh();
    
    /*!
     * Mesh scaling for solver (gamma factor)
     * (in particle rest frame, the longitudinal length enlarged)
     */
    const Vector_t& getMeshScaling() const;
    
    Vektor<int, 3> getBaseLevelGridPoints() const;
    
    const int& maxLevel() const;
    const int& finestLevel() const;
    
    /*!
     * @returns the time of the bunch [s]
     */
    double getT() const;
    
    void redistributeGrids(int how);
    
    
protected:
    /*
     * AmrMesh functions
     */
    
    /*!
     * Update the grids and the distributionmapping for a specific level
     * (inherited from AmrMesh)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void RemakeLevel (int lev, AmrReal_t time,
                      const AmrGrid_t& new_grids, const AmrProcMap_t& new_dmap);
    
    /*!
     * Create completeley new grids for a level (inherited from AmrMesh)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void MakeNewLevel (int lev, AmrReal_t time,
                       const AmrGrid_t& new_grids, const AmrProcMap_t& new_dmap);
    
    /*!
     * Clean up a level.
     * @param lev to free allocated memory.
     */
    void ClearLevel(int lev);
    
    /*!
     * Is called in the AmrMesh function for performing tagging. (inherited from AmrMesh)
     */
    virtual void ErrorEst(int lev, TagBoxArray_t& tags,
                          AmrReal_t time, int ngrow) override;
    
    
    /*!
     * Make a new level from scratch using provided BoxArray and
     * DistributionMapping.
     * Only used during initialization.
     * 
     * Remark: Not used in OPAL
     * 
     * @param lev to create
     * @param time of simulation
     * @param ba the boxes
     * @param dm the grid distribution among cores
     */
    void MakeNewLevelFromScratch (int lev, AmrReal_t time,
                                  const AmrGrid_t& ba,
                                  const AmrProcMap_t& dm);

    /*!
     * Make a new level using provided BoxArray and
     * DistributionMapping and fill with interpolated coarse level data.
     * 
     * Remark: Not used in OPAL
     * 
     * @param lev to create
     * @param time of simulation
     * @param ba the boxes
     * @param dm the grid distribution among cores
     */
    void MakeNewLevelFromCoarse (int lev, AmrReal_t time,
                                 const AmrGrid_t& ba,
                                 const AmrProcMap_t& dm);

private:
    /*!
     * Update of mesh according to chosen refinement strategy.
     * Rebuild levels finer than lbase.
     * @param lbase base level to start regriding
     * @param time of regrid
     */
    void doRegrid_m(int lbase, double time);
    
    /*!
     * Called within doRegrid_m(). Redistribute
     * particles and delete levels on particle side.
     * @param old_finest level
     */
    void postRegrid_m(int old_finest);

    
    /* ATTENTION
     * The tagging routines assume the particles to be in the
     * AMR domain, i.e. [-1, 1]^3
     */
    
    /*!
     * Mark a cell for refinement if the value is greater equal
     * than some amount of charge (AmrObject::chargedensity_m).
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForChargeDensity_m(int lev, TagBoxArray_t& tags,
                               AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if the potential value is greater
     * equal than the maximum value of the potential on the grid
     * scaled by some factor [0, 1] (AmrObject::scaling_m)
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForPotentialStrength_m(int lev, TagBoxArray_t& tags,
                                   AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if one of the electric field components
     * is greater equal the maximum electric field value per direction scaled
     * by some factor [0, 1] (AmrObject::scaling_m).
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForEfield_m(int lev, TagBoxArray_t& tags,
                        AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if at least one particle has
     * a high momentum. The lower bound is specified by the
     * maximum momenta per level scaled by some factor [0, 1]
     * (AmrObject::scaling_m).
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForMomenta_m(int lev, TagBoxArray_t& tags,
                         AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if it contains at most
     * AmrObject::maxNumPart_m particles.
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForMaxNumParticles_m(int lev, TagBoxArray_t& tags,
                                 AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if it contains at least
     * AmrObject::minNumPart_m particles.
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForMinNumParticles_m(int lev, TagBoxArray_t& tags,
                                 AmrReal_t time, int ngrow);
    
    /*!
     * Use particle BoxArray and DistributionMapping for AmrObject and
     * reset geometry for bunch
     * 
     * @param nGridPts per dimension (nx, ny, nz / nt)
     */
    void initBaseLevel_m(const AmrIntArray_t& nGridPts);
    
    /*!
     * AMReX uses the ParmParse object to initialize
     * parameters like the maximum level etc.
     * This function initializes "all" of them.
     * 
     * @param info all parameters that we set over the
     *             OPAL input file
     * @param layout_p of bunch
     */
    static void initParmParse_m(const AmrInfo& info, AmrLayout_t* layout_p);
    
    /*!
     * Fill the physical / mesh boundary values
     * 
     * @param mf field to fill boundary values
     */
    void fillPhysbc_m(AmrField_t& mf, int lev = 0);
    
    
//     void gradient(int lev) {
//         
//         phi_m[lev]->FillBoundary(geom[lev].periodicity());
//         
//         for (MFIter_t mfi(*phi_m[lev], true); mfi.isValid(); ++mfi) {
//             const Box_t&  tilebx  = mfi.tilebox();
//             FArrayBox_t&  fab = (*phi_m[lev])[mfi];
//             FArrayBox_t& efab = (*efield_m[lev])[mfi];
//             
//             const int*  tlo     = tilebx.loVect();
//             const int*  thi     = tilebx.hiVect();
//             
//             for (int i = tlo[0]; i <= thi[0]; ++i) {
//                 for (int j = tlo[1]; j <= thi[1]; ++j) {
//                     for (int k = tlo[2]; k <= thi[2]; ++k) {
//                         
//                         amrex::IntVect iv(D_DECL(i,j,k));
//                         
//                         // x-field
//                         amrex::IntVect liv(D_DECL(i-1,j,k));
//                         amrex::IntVect riv(D_DECL(i+1,j,k));
//                         efab(iv, 0) = 0.5 * (fab(liv) - fab(riv));
//                         
//                         // y-field
//                         liv = amrex::IntVect(D_DECL(i,j-1,k));
//                         riv = amrex::IntVect(D_DECL(i,j+1,k));
//                         efab(iv, 1) = 0.5 * (fab(liv) - fab(riv));
//                         
//                         // z-field
//                         liv = amrex::IntVect(D_DECL(i,j-1,k));
//                         riv = amrex::IntVect(D_DECL(i,j+1,k));
//                         efab(iv, 2) = 0.5 * (fab(liv) - fab(riv));
//                     }
//                 }
//             }
//         }
//     }
    
private:
    /// bunch used for tagging strategies
    AmrPartBunch *bunch_mp;
    
    // the layout of the bunch
    AmrLayout_t  *layout_mp;
    
    /// charge density on the grid for all levels
    AmrScalarFieldContainer_t rho_m;
    
    /// scalar potential on the grid for all levels
    AmrScalarFieldContainer_t phi_m;
    
    /// vector field on the grid for all levels
    AmrVectorFieldContainer_t efield_m;
    
    /// in particle rest frame, the longitudinal length enlarged
    Vector_t meshScaling_m;

    // only used in case of potential and efield tagging
    std::vector<bool> isFirstTagging_m;

    bool isPoissonSolved_m;
};

#endif
