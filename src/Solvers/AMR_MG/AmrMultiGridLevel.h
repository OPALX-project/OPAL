#ifndef AMR_MULTI_GRID_LEVEL
#define AMR_MULTI_GRID_LEVEL

#include <vector>

#include <AMReX_IntVect.H>

#include "Ippl.h"

#include "AmrMultiGridDefs.h"

#include <unordered_map>

template <class MatrixType, class VectorType>
class AmrMultiGridLevel {
    
public:
    typedef amr::AmrField_t             AmrField_t;
    typedef amr::AmrGeometry_t          AmrGeometry_t;
    typedef std::unique_ptr<AmrField_t> AmrField_u;
    typedef std::shared_ptr<AmrField_t> AmrField_s;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    typedef MatrixType                  matrix_t;
    typedef VectorType                  vector_t;
    typedef amrex::BaseFab<int>         basefab_t;
    typedef amrex::FabArray<basefab_t>  mask_t;
    typedef std::shared_ptr<AmrBoundary<AmrMultiGridLevel<MatrixType,
                                                          VectorType
                                                          >
                                        >
                            > boundary_t;
    
    typedef amr::comm_t comm_t;
    typedef amr::dmap_t dmap_t;
    typedef amr::node_t node_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::local_ordinal_t lo_t;
    
    /// Type for matrix indices
    typedef std::vector<go_t> indices_t;
    
    /// Type for matrix entries
    typedef std::vector<scalar_t> coefficients_t;
    
    // Type with matrix index (column) and coefficient value
    typedef std::unordered_map<go_t, scalar_t> umap_t;
    
    // covered   : ghost cells covered by valid cells of this FabArray
    //             (including periodically shifted valid cells)
    // notcovered: ghost cells not covered by valid cells
    //             (including ghost cells outside periodic boundaries) (BNDRY)
    // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
    // interior  : interior cells (i.e., valid cells)
    enum Mask {
        COVERED   = -1,
        INTERIOR  =  0,
        BNDRY     =  1,
        PHYSBNDRY =  2
    };

    // NO   : not a refined cell 
    // YES  : cell got refined
    enum Refined {
        YES = 0,
        NO  = 1
    };
    
public:
    /*!
     * @param mesh scaling due to particle rest frame
     * @param _grids of this level
     * @param _dmap AMReX core distribution map
     * @param _geom of domain
     * @param rr refinement ratio
     * @param bc physical boundaries (x, y, z)
     * @param comm MPI communicator
     * @param node Kokkos node type (Serial, OpenMP, CUDA)
     */
    AmrMultiGridLevel(const Vector_t& meshScaling,
                      const amrex::BoxArray& _grids,
                      const amrex::DistributionMapping& _dmap,
                      const AmrGeometry_t& _geom,
                      const AmrIntVect_t& rr,
                      const boundary_t* bc,
                      const Teuchos::RCP<comm_t>& comm,
                      const Teuchos::RCP<node_t>& node);
    
    ~AmrMultiGridLevel();
    
    /*!
     * Map a 2D / 3D grid point to an array index
     * @param iv grid point (i, j, k)
     */
    go_t serialize(const AmrIntVect_t& iv) const;
    
    /*!
     * Checks if grid point is on the physical / mesh boundary
     * @param iv grid point (i, j, k)
     */
    bool isBoundary(const AmrIntVect_t& iv) const;
    
    /*!
     * Checks all directions if physical / mesh boundary.
     * @param iv is the cell where we want to have the boundary value
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     */
    bool applyBoundary(const AmrIntVect_t& iv,
                       umap_t& map,
                       const scalar_t& value);
    
    /*!
     * Slightly faster version of apply().
     * @param iv is the cell where we want to have the boundary value
     * @param fab is the mask
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     * @precondition Basefab needs to be a mask with
     * AmrMultiGridLevel::Mask::PHYSBNDRY
     */
    bool applyBoundary(const AmrIntVect_t& iv,
                       const basefab_t& fab,
                       umap_t& map,
                       const scalar_t& value);
    
    /*!
     * Apply boundary in a certain direction.
     * @param iv is the cell where we want to have the boundary value
     * @param dir direction of physical / mesh boundary
     * @param map with indices global matrix indices and matrix values
     * @param value matrix entry (coefficients)
     */
    void applyBoundary(const AmrIntVect_t& iv,
                       const lo_t& dir,
                       umap_t& map,
                       const scalar_t& value);
    
    const AmrIntVect_t& refinement() const;
    
    /*!
     * @returns the mesh spacing in particle rest frame
     */
    const scalar_t* cellSize() const;
    
    /*!
     * @returns the mesh spacing in particle rest frame for a
     * certain direction
     */
    const scalar_t& cellSize(lo_t dir) const;
    
    /*!
     * @returns the inverse mesh spacing in particle rest frame
     */
    const scalar_t* invCellSize() const;
    
    /*!
     * @returns the inverse mesh spacing in particle rest
     * frame for a certain direction
     */
    const scalar_t& invCellSize(lo_t dir) const;

    /*!
     * Check if a cell on that level is valid
     * @param iv is the cell
     * @returns true if the cell is valid
     */
    bool isValid(const AmrIntVect_t& iv) const;

private:
    /*!
     * Build a mask specifying if a grid point is covered,
     * an interior cell, at physical boundary or at interior boundary
     */
    void buildLevelMask_m();
    
    /*!
     * Build Tpetra::Map of this level
     * @param comm MPI communicator
     * @param node Kokkos node type
     */
    void buildMap_m(const Teuchos::RCP<comm_t>& comm,
                    const Teuchos::RCP<node_t>& node);
    
public:
    const amrex::BoxArray& grids;           ///< boxes of this level
    const amrex::DistributionMapping& dmap; ///< AMReX core distribution map
    const AmrGeometry_t& geom;              ///< geometry of this problem
    
    Teuchos::RCP<dmap_t> map_p;         ///< Tpetra core map
    
    Teuchos::RCP<matrix_t> Anf_p;       ///< no fine Poisson matrix
    Teuchos::RCP<matrix_t> R_p;         ///< restriction matrix
    Teuchos::RCP<matrix_t> I_p;         ///< interpolation matrix
    Teuchos::RCP<matrix_t> Bcrse_p;     ///< boundary from coarse cells
    Teuchos::RCP<matrix_t> Bfine_p;     ///< boundary from fine cells
    Teuchos::RCP<matrix_t> Awf_p;       ///< composite Poisson matrix
    
    /// gradient matrices in x, y, and z to compute electric field
    Teuchos::RCP<matrix_t> G_p[AMREX_SPACEDIM];
    
    Teuchos::RCP<vector_t> rho_p;       ///< charge density
    Teuchos::RCP<vector_t> phi_p;       ///< potential vector
    Teuchos::RCP<vector_t> residual_p;  ///< residual over all cells
    Teuchos::RCP<vector_t> error_p;     ///< error over all cells
    Teuchos::RCP<matrix_t> UnCovered_p; ///< uncovered cells
    
    std::unique_ptr<mask_t> mask;       ///< interior, phys boundary, interface, covered
    std::unique_ptr<mask_t> refmask;    ///< covered (i.e. refined) or not-covered
    std::unique_ptr<mask_t> crsemask;
    
private:
    go_t nr_m[AMREX_SPACEDIM];          ///< number of grid points
    
    AmrIntVect_t rr_m;                  ///< refinement
    
    boundary_t bc_mp[AMREX_SPACEDIM];   ///< boundary conditions
    
    scalar_t dx_m[AMREX_SPACEDIM];      ///< cell size in particle rest frame
    scalar_t invdx_m[AMREX_SPACEDIM];   ///< inverse cell size in particle rest frame
};


#include "AmrMultiGridLevel.hpp"

#endif
