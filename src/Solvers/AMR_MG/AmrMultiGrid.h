#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "AmrMultiGridCore.h"

#include "Solvers/AmrPoissonSolver.h"
#include "Amr/AmrBoxLib.h"

#include "AmrMultiGridLevel.h"

#include <fstream>

#define AMR_MG_TIMER true
#define AMR_MG_WRITE false

class AmrMultiGrid : public AmrPoissonSolver< AmrBoxLib > {
    
public:
    typedef amr::matrix_t         matrix_t;
    typedef amr::vector_t         vector_t;
    typedef amr::multivector_t    mv_t;
    typedef amr::dmap_t           dmap_t;
    typedef amr::comm_t           comm_t;
    typedef amr::local_ordinal_t  lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::scalar_t         scalar_t;
    
    typedef AmrMultiGridLevel<matrix_t, vector_t> AmrMultiGridLevel_t;
    
    typedef AmrMultiGridLevel_t::AmrField_t     AmrField_t;
    typedef AmrMultiGridLevel_t::AmrGeometry_t  AmrGeometry_t;
    typedef AmrMultiGridLevel_t::AmrField_u     AmrField_u;
    typedef AmrMultiGridLevel_t::AmrField_s     AmrField_s;
    typedef AmrMultiGridLevel_t::AmrIntVect_t   AmrIntVect_t;
    typedef AmrMultiGridLevel_t::indices_t      indices_t;
    typedef AmrMultiGridLevel_t::coefficients_t coefficients_t;
    typedef AmrMultiGridLevel_t::umap_t         umap_t;
    typedef AmrMultiGridLevel_t::boundary_t     boundary_t;
    
    typedef BottomSolver<
        Teuchos::RCP<matrix_t>,
        Teuchos::RCP<mv_t>,
        AmrMultiGridLevel_t
    > bsolver_t;
    
    typedef BelosBottomSolver<AmrMultiGridLevel_t>      BelosSolver_t;
    typedef Amesos2BottomSolver<AmrMultiGridLevel_t>    Amesos2Solver_t;
    typedef MueLuBottomSolver<AmrMultiGridLevel_t>      MueLuSolver_t;
    
    typedef AmrPreconditioner<matrix_t, AmrMultiGridLevel_t> preconditioner_t;
    
    typedef Ifpack2Preconditioner<AmrMultiGridLevel_t> Ifpack2Preconditioner_t;
    typedef MueLuPreconditioner<AmrMultiGridLevel_t> MueLuPreconditioner_t;
    
    typedef amrex::BoxArray boxarray_t;
    typedef amrex::Box box_t;
    typedef amrex::BaseFab<int> basefab_t;
    typedef amrex::FArrayBox farraybox_t;
            
    typedef AmrSmoother::Smoother Smoother;
    
    typedef amr::Preconditioner Preconditioner;
    
    /// Supported interpolaters for prolongation operation
    enum Interpolater {
        TRILINEAR = 0,
        LAGRANGE,
        PIECEWISE_CONST
    };
    
    /// Supported bottom solvers
    enum BaseSolver {
        // all Belos
        BICGSTAB,
        MINRES,
        PCPG,
        CG,
        GMRES,
        STOCHASTIC_CG,
        RECYCLING_CG,
        RECYCLING_GMRES
        // all Amesos2
#ifdef HAVE_AMESOS2_KLU2
        , KLU2
#endif
#ifdef HAVE_AMESOS2_SUPERLU
        , SUPERLU
#endif
#ifdef HAVE_AMESOS2_UMFPACK  
        , UMFPACK
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
        , PARDISO_MKL
#endif
#ifdef HAVE_AMESOS2_MUMPS
        , MUMPS
#endif
#ifdef HAVE_AMESOS2_LAPACK
        , LAPACK
#endif
        // all MueLu
        , SA
        // add others ...
    };
    
    /// Supported physical boundaries
    enum Boundary {
        DIRICHLET = 0,
        OPEN,
        PERIODIC
    };
    
    /// Supported convergence criteria
    enum Norm {
        L1,
        L2,
        LINF
    };
    
public:
    
    /*!
     * Instantiation used in Structure/FieldSolver.cpp
     * @param itsAmrObject_p has information about refinemen ratios, etc.
     * @param bsolver bottom solver
     * @param prec preconditioner for bottom solver
     * @param rebalance of preconditioner (SA only)
     * @param bcx boundary condition in x
     * @param bcy boundary condition in y
     * @param bcz boundary condition in z
     * @param smoother for level solution
     * @param nSweeps when smoothing
     * @param interp interpolater between levels
     * @param norm for convergence criteria
     */
    AmrMultiGrid(AmrBoxLib* itsAmrObject_p,
                 const std::string& bsolver,
                 const std::string& prec,
                 const bool& rebalance,
                 const std::string& reuse,
                 const std::string& bcx,
                 const std::string& bcy,
                 const std::string& bcz,
                 const std::string& smoother,
                 const std::size_t& nSweeps,
                 const std::string& interp,
                 const std::string& norm);
    
    /*!
     * Used in OPAL
     * 
     * @param rho right-hand side charge density on grid [C / m]
     * @param phi electrostatic potential (unknown) [V]
     * @param efield electric field [V / m]
     * @param baseLevel for solve
     * @param finestLevel for solve
     * @param prevAsGuess use of previous solution as initial guess
     */
    void solve(AmrScalarFieldContainer_t &rho,
               AmrScalarFieldContainer_t &phi,
               AmrVectorFieldContainer_t &efield,
               unsigned short baseLevel,
               unsigned short finestLevel,
               bool prevAsGuess = true);
    
    /*!
     * Specify the number of smoothing steps
     * @param nSweeps for each smoothing step
     */
    void setNumberOfSweeps(const std::size_t& nSweeps);
    
    /*!
     * Specify the maximum number of iterations
     * @param maxiter \f$ [0, \infty[ \f$
     */
    void setMaxNumberOfIterations(const std::size_t& maxiter);
    
    /*!
     * Obtain some convergence info
     * @returns the number of iterations till convergence
     */
    std::size_t getNumIters();
    
    /*!
     * Obtain the residual norm of a level
     * @param level for which error is requested
     * @returns the norm of the residual
     */
    scalar_t getLevelResidualNorm(lo_t level);
    
    /*!
     * Enable solver info dumping into SDDS file
     */
    void setVerbose(bool verbose);
    
    double getXRangeMin(unsigned short level = 0);
    double getXRangeMax(unsigned short level = 0);
    double getYRangeMin(unsigned short level = 0);
    double getYRangeMax(unsigned short level = 0);
    double getZRangeMin(unsigned short level = 0);
    double getZRangeMax(unsigned short level = 0);
    
    /**
     * Print information abour tolerances.
     * @param os output stream where to write to
     */
    Inform &print(Inform &os) const;
    
private:
    
    /*!
     * Instantiate boundary object
     * @param bc boundary conditions
     * @precondition length must be equal to AMREX_SPACEDIM
     */
    void initPhysicalBoundary_m(const Boundary* bc);
    
    /*!
     * Instantiate all levels and set boundary conditions
     * @param rho is the charge density
     * @param geom is the geometry
     * @param regrid was performed
     */
    void initLevels_m(const amrex::Vector<AmrField_u>& rho,
                      const amrex::Vector<AmrGeometry_t>& geom,
                      bool regrid);
    
    /*!
     * Clear masks (required to build matrices) no longer needed.
     */
    void clearMasks_m();
    
    /*!
     * Reset potential to zero (currently)
     * @param reset solution to initial guess (zero)
     */
    void initGuess_m(bool reset);
    
    /*!
     * Actual solve.
     * @returns the the max. residual
     */
    scalar_t iterate_m();
    
    /*!
     * Compute norms / level and check convergence
     * @returns true if converged
     */
    bool isConverged_m(std::vector<scalar_t>& rhsNorms,
                       std::vector<scalar_t>& resNorms);
    
    /*!
     * Compute composite residual of a level
     * @param r is the residual to compute
     * @param b is the right-hand side
     * @param x is the left-hand side
     * @param level to solve for
     */
    void residual_m(const lo_t& level,
                    Teuchos::RCP<vector_t>& r,
                    const Teuchos::RCP<vector_t>& b,
                    const Teuchos::RCP<vector_t>& x);
    
    /*!
     * Recursive call.
     * @param level to relax
     */
    void relax_m(const lo_t& level);
    
    /*!
     * Compute the residual of a level without considering refined level.
     * @param result is computed
     * @param rhs is the right-hand side
     * @param crs_rhs is the coarse right-hand side for internal boundary
     * @param b is the left-hand side
     * @param level to solver for
     */
    void residual_no_fine_m(const lo_t& level,
                            Teuchos::RCP<vector_t>& result,
                            const Teuchos::RCP<vector_t>& rhs,
                            const Teuchos::RCP<vector_t>& crs_rhs,
                            const Teuchos::RCP<vector_t>& b);
    
#if AMR_MG_WRITE
    /*!
     * Dumps the residual norm per level into a file (for each iteration).
     */
    void writeResidualNorm_m();
#endif
    
    /*!
     * Vector norm computation.
     * @param x is the vector for which we compute the norm
     * @returns the evaluated norm of a level
     */
    scalar_t evalNorm_m(const Teuchos::RCP<const vector_t>& x);
    
    /*!
     * Initial convergence criteria values.
     * @param rhsNorms per level of right-hand side (is filled)
     * @param resNorms per level of residual (is filled)
     */
    void initResidual_m(std::vector<scalar_t>& rhsNorms,
                        std::vector<scalar_t>& resNorms);
    
    /*!
     * @param efield to compute
     */
    void computeEfield_m(AmrVectorFieldContainer_t& efield);
    
    /*!
     * Build all matrices and vectors, i.e. AMReX to Trilinos
     * @param matrices if we need to build matrices as well or only vectors.
     */
    void setup_m(const amrex::Vector<AmrField_u>& rho,
                 const amrex::Vector<AmrField_u>& phi,
                 const bool& matrices = true);
    
    /*!
     * Build all matrices and vectors needed for single-level computation
     */
    void buildSingleLevel_m(const amrex::Vector<AmrField_u>& rho,
                            const amrex::Vector<AmrField_u>& phi,
                            const bool& matrices = true);
    
    /*!
     * Build all matrices and vectors needed for multi-level computation
     */
    void buildMultiLevel_m(const amrex::Vector<AmrField_u>& rho,
                           const amrex::Vector<AmrField_u>& phi,
                           const bool& matrices = true);

    /*!
     * Set matrix and vector pointer
     * @param level for which we fill matrix + vector
     * @param matrices if we need to set matrices as well or only vectors.
     */
    void open_m(const lo_t& level, const bool& matrices);
    
    /*!
     * Call fill complete
     * @param level for which we filled matrix
     * @param matrices if we set matrices as well.
     */
    void close_m(const lo_t& level, const bool& matrices);
    
    /*!
     * Build the Poisson matrix for a level assuming no finer level (i.e. the whole fine mesh
     * is taken into account).
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param level for which we build the Poisson matrix
     */
    void buildNoFinePoissonMatrix_m(const lo_t& level,
                                    const go_t& gidx,
                                    const AmrIntVect_t& iv,
                                    const basefab_t& mfab,
                                    const scalar_t* invdx2);
    
    /*!
     * Build the Poisson matrix for a level that got refined (it does not take the covered
     * cells (covered by fine cells) into account). The finest level does not build such a matrix.
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param rfab is the mask between levels
     * @param level for which we build the special Poisson matrix
     */
    void buildCompositePoissonMatrix_m(const lo_t& level,
                                       const go_t& gidx,
                                       const AmrIntVect_t& iv,
                                       const basefab_t& mfab,
                                       const basefab_t& rfab,
                                       const basefab_t& cfab,
                                       const scalar_t* invdx2);
    
    /*!
     * Build a matrix that averages down the data of the fine cells down to the
     * corresponding coarse cells. The base level does not build such a matrix.
     * \f[
     *      x^{(l)} = R\cdot x^{(l+1)}
     * \f]
     * @param iv is the current cell
     * @param rfab is the mask between levels
     * @param level for which to build restriction matrix
     */
    void buildRestrictionMatrix_m(const lo_t& level,
                                  const go_t& gidx,
                                  const AmrIntVect_t& iv,
                                  D_DECL(const go_t& ii,
                                         const go_t& jj,
                                         const go_t& kk),
                                  const basefab_t& rfab);
    
    /*!
     * Interpolate data from coarse cells to appropriate refined cells. The interpolation
     * scheme is allowed only to have local cells in the stencil, i.e.
     * 2D --> 4, 3D --> 8.
     * \f[
     *      x^{(l)} = I\cdot x^{(l-1)}
     * \f]
     * @param iv is the current cell
     * @param level for which to build the interpolation matrix. The finest level
     * does not build such a matrix.
     */
    void buildInterpolationMatrix_m(const lo_t& level,
                                    const go_t& gidx,
                                    const AmrIntVect_t& iv,
                                    const basefab_t& cfab);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the fine boundary values from the coarse values, i.e.
     * \f[
     *      x^{(l)} = B_{crse}\cdot x^{(l-1)}
     * \f]
     * Dirichlet boundary condition
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param cells all fine cells that are at the crse-fine interface
     * @param level the base level is omitted
     */
    void buildCrseBoundaryMatrix_m(const lo_t& level,
                                   const go_t& gidx,
                                   const AmrIntVect_t& iv,
                                   const basefab_t& mfab,
                                   const basefab_t& cfab,
                                   const scalar_t* invdx2);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the coarse boundary values from fine values, i.e.
     * \f[
     *      x^{(l)} = B_{fine}\cdot x^{(l+1)}
     * \f]
     * Dirichlet boundary condition. Flux matching.
     * @param iv is the current cell
     * @param cells all coarse cells that are at the crse-fine interface but are
     * not refined
     * @param crse_fine_ba coarse cells that got refined
     * @param level the finest level is omitted
     */
    void buildFineBoundaryMatrix_m(const lo_t& level,
                                   const go_t& gidx,
                                   const AmrIntVect_t& iv,
                                   const basefab_t& mfab,
                                   const basefab_t& rfab,
                                   const basefab_t& cfab);
    
    /*!
     * Copy data from AMReX to Trilinos
     * @param rho is the charge density
     * @param level for which to copy
     */
    void buildDensityVector_m(const lo_t& level,
                              const AmrField_t& rho);
    
    /*!
     * Copy data from AMReX to Trilinos
     * @param phi is the potential
     * @param level for which to copy
     */
    void buildPotentialVector_m(const lo_t& level,
                                const AmrField_t& phi);
    
    /*!
     * Gradient matrix is used to compute the electric field
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param level for which to compute
     */
    void buildGradientMatrix_m(const lo_t& level,
                               const go_t& gidx,
                               const AmrIntVect_t& iv,
                               const basefab_t& mfab,
                               const scalar_t* invdx);
    
    /*!
     * Data transfer from AMReX to Trilinos.
     * @param mf is the multifab of a level
     * @param comp component to copy
     * @param mv is the vector to be filled
     * @param level where to perform
     */
    void amrex2trilinos_m(const lo_t& level,
                          const lo_t& comp,
                          const AmrField_t& mf,
                          Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Data transfer from Trilinos to AMReX.
     * @param level to copy
     * @param comp component to copy
     * @param mf is the multifab to be filled
     * @param mv is the corresponding Trilinos vector
     */
    void trilinos2amrex_m(const lo_t& level,
                          const lo_t& comp,
                          AmrField_t& mf,
                          const Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Some indices might occur several times due to boundary conditions, etc. We
     * avoid this by filling a map and then copying the data to a vector for filling
     * the matrices. The map gets cleared inside the function.
     * @param indices in matrix
     * @param values are the coefficients
     */
    void map2vector_m(umap_t& map, indices_t& indices,
                      coefficients_t& values);
    
    /*!
     * Perform one smoothing step
     * @param e error to update (left-hand side)
     * @param r residual (right-hand side)
     * @param level on which to relax
     */
    void smooth_m(const lo_t& level,
                  Teuchos::RCP<vector_t>& e,
                  Teuchos::RCP<vector_t>& r);
    
    /*!
     * Restrict coarse level residual based on fine level residual
     * @param level to restrict
     */
    void restrict_m(const lo_t& level);
    
    /*!
     * Update error of fine level based on error of coarse level
     * @param level to update
     */
    void prolongate_m(const lo_t& level);
    
    /*!
     * Average data from fine level to coarse
     * @param level finest level is omitted
     */
    void averageDown_m(const lo_t& level);
    
    /*!
     * Instantiate interpolater
     * @param interp interpolater type
     */
    void initInterpolater_m(const Interpolater& interp);
    
    /*!
     * Instantiate interface interpolater
     * @param interface handler
     */
    void initCrseFineInterp_m(const Interpolater& interface);
    
    /*!
     * Instantiate a bottom solver
     * @param solver type
     * @param rebalance solver (SA only)
     * @param reuse types of SA hierarchy
     */
    void initBaseSolver_m(const BaseSolver& solver,
                          const bool& rebalance,
                          const std::string& reuse);
    
    /*!
     * Instantiate a preconditioner for the bottom solver
     * @param precond type
     * @param rebalance preconditioner (SA only)
     * @param reuse types of SA hierarchy
     */
    void initPrec_m(const Preconditioner& prec,
                    const bool& rebalance,
                    const std::string& reuse);
    
    /*!
     * Convertstring to enum Boundary
     * @param bc boundary condition
     */
    Boundary convertToEnumBoundary_m(const std::string& bc);
    
    /*!
     * Converts string to enum Interpolater
     * @param interp interpolater
     */
    Interpolater convertToEnumInterpolater_m(const std::string& interp);
    
    /*!
     * Converts string to enum BaseSolver
     * @param bsolver bottom solver
     */
    BaseSolver convertToEnumBaseSolver_m(const std::string& bsolver);
    
    /*!
     * Converts string to enum Preconditioner
     * @param prec preconditioner
     */
    Preconditioner convertToEnumPreconditioner_m(const std::string& prec);
    
    /*!
     * Converts string to enum Smoother
     * @param smoother of level solution
     */
    Smoother convertToEnumSmoother_m(const std::string& smoother);
    
    /*!
     * Converts string to enum Norm
     * @param norm either L1, L2, LInf
     */
    Norm convertToEnumNorm_m(const std::string& norm);
    
    /*!
     * SDDS header is written by root core
     * @param outfile output stream
     */
    void writeSDDSHeader_m(std::ofstream& outfile);
    
    /*!
     * SDDS data write (done by root core)
     * @param error to write
     */
    void writeSDDSData_m(const scalar_t& error);
    
#if AMR_MG_TIMER
    /*!
     * Create timers
     */
    void initTimer_m();
#endif
    
private:
    Teuchos::RCP<comm_t> comm_mp;       ///< communicator
    Teuchos::RCP<amr::node_t> node_mp;  ///< kokkos node
    
    /// interpolater without coarse-fine interface
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interp_mp;
    
    /// interpolater for coarse-fine interface
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interface_mp;
    
    std::size_t nIter_m;            ///< number of iterations till convergence
    std::size_t bIter_m;            ///< number of iterations of bottom solver
    std::size_t maxiter_m;          ///< maximum number of iterations allowed
    std::size_t nSweeps_m;          ///< number of smoothing iterations
    Smoother smootherType_m;        ///< type of smoother
    
    /// container for levels
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
    /// bottom solver
    std::shared_ptr<bsolver_t> solver_mp;
    
    /// error smoother
    std::vector<std::shared_ptr<AmrSmoother> > smoother_m;
    
    /// preconditioner for bottom solver
    std::shared_ptr<preconditioner_t> prec_mp;
    
    int lbase_m;            ///< base level (currently only 0 supported)
    int lfine_m;            ///< fineste level
    int nlevel_m;           ///< number of levelss
    
    boundary_t bc_m[AMREX_SPACEDIM];    ///< boundary conditions
    int nBcPoints_m;                    ///< maximum number of stencils points for BC
    
    Norm norm_m;            ///< norm for convergence criteria (l1, l2, linf)
    std::string snorm_m;    ///< norm for convergence criteria
    
    const scalar_t eps_m;   ///< rhs scale for convergence
    
    bool verbose_m;                 ///< If true, a SDDS file is written
    std::string fname_m;            ///< SDDS filename
    std::ios_base::openmode flag_m; ///< std::ios::out or std::ios::app
    
#if AMR_MG_TIMER
    IpplTimings::TimerRef buildTimer_m;         ///< timer for matrix and vector construction
    IpplTimings::TimerRef restrictTimer_m;      ///< timer for restriction operation
    IpplTimings::TimerRef smoothTimer_m;        ///< timer for all smoothing steps
    IpplTimings::TimerRef interpTimer_m;        ///< prolongation timer
    IpplTimings::TimerRef efieldTimer_m;        ///< timer for e-field calculation + copy
    IpplTimings::TimerRef averageTimer_m;       ///> timer for average down process
    IpplTimings::TimerRef bottomTimer_m;        ///< bottom solver timer
    IpplTimings::TimerRef dumpTimer_m;          ///< write SDDS file timer
#endif
};


inline Inform &operator<<(Inform &os, const AmrMultiGrid &fs) {
    return fs.print(os);
}

#endif
