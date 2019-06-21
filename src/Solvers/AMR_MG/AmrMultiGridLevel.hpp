#define AMR_NO_SCALE false


template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType,
                  VectorType>::AmrMultiGridLevel(const Vector_t& meshScaling,
                                                 const amrex::BoxArray& _grids,
                                                 const amrex::DistributionMapping& _dmap,
                                                 const AmrGeometry_t& _geom,
                                                 const AmrIntVect_t& rr,
                                                 const boundary_t* bc,
                                                 const Teuchos::RCP<comm_t>& comm,
                                                 const Teuchos::RCP<node_t>& node)
    : grids(_grids),
      dmap(_dmap),
      geom(_geom),
      map_p(Teuchos::null),
      Anf_p(Teuchos::null),
      R_p(Teuchos::null),
      I_p(Teuchos::null),
      Bcrse_p(Teuchos::null),
      Bfine_p(Teuchos::null),
      Awf_p(Teuchos::null),
      rho_p(Teuchos::null),
      phi_p(Teuchos::null),
      residual_p(Teuchos::null),
      error_p(Teuchos::null),
      UnCovered_p(Teuchos::null),
      refmask(nullptr),
      crsemask(nullptr),
      rr_m(rr)
{
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        G_p[j] = Teuchos::null;
        
        nr_m[j] = _geom.Domain().length(j);
        
#if AMR_NO_SCALE
        // mesh spacing in particle rest frame
        dx_m[j] = geom.CellSize(j);
        invdx_m[j] = geom.InvCellSize(j);
#else
        // mesh spacing in particle rest frame
        dx_m[j] = meshScaling[j] * geom.CellSize(j);
        invdx_m[j] = meshScaling[j] * geom.InvCellSize(j);
#endif
        
        bc_mp[j] = bc[j];
    }
    
    this->buildLevelMask_m();
    
    this->buildMap_m(comm, node);
    
    
    residual_p = Teuchos::rcp( new vector_t(map_p, false) );
    error_p = Teuchos::rcp( new vector_t(map_p, false) );
}


template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType, VectorType>::~AmrMultiGridLevel()
{
    map_p = Teuchos::null;
    
    Anf_p = Teuchos::null;
    R_p = Teuchos::null;
    I_p = Teuchos::null;
    Bcrse_p = Teuchos::null;
    Bfine_p = Teuchos::null;
    Awf_p = Teuchos::null;
    
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        G_p[j] = Teuchos::null;
    
    UnCovered_p = Teuchos::null;
    
    rho_p = Teuchos::null;
    phi_p = Teuchos::null;
    residual_p = Teuchos::null;
    error_p = Teuchos::null;
}


template <class MatrixType, class VectorType>
typename AmrMultiGridLevel<MatrixType, VectorType>::go_t
AmrMultiGridLevel<MatrixType, VectorType>::serialize(const AmrIntVect_t& iv) const {
#if AMREX_SPACEDIM == 3
    return iv[0] + (iv[1] + nr_m[1] * iv[2]) * nr_m[0];
#else
    return iv[0] + iv[1] * nr_m[0];
#endif
}


template <class MatrixType, class VectorType>
bool AmrMultiGridLevel<MatrixType, VectorType>::isBoundary(const AmrIntVect_t& iv) const {
    // it doesn't matter with which direction we check, since it checks all
    return bc_mp[0]->isBoundary(iv, &nr_m[0]);
}


template <class MatrixType, class VectorType>
bool AmrMultiGridLevel<MatrixType, VectorType>::applyBoundary(const AmrIntVect_t& iv,
                                                              umap_t& map,
                                                              const scalar_t& value)
{
    bool applied = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ( bc_mp[d]->isBoundary(iv, d, &nr_m[0]) ) {
            applied = true;
            bc_mp[d]->apply(iv, d, map, value, this, &nr_m[0]);
        }
    }
    return applied;
}


template <class MatrixType, class VectorType>
bool AmrMultiGridLevel<MatrixType, VectorType>::applyBoundary(const AmrIntVect_t& iv,
                                                              const basefab_t& fab,
                                                              umap_t& map,
                                                              const scalar_t& value)
{
    if ( fab(iv) != Mask::PHYSBNDRY )
        return false;
    
    bool applied = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ( bc_mp[d]->isBoundary(iv, d, &nr_m[0]) ) {
            applied = true;
            bc_mp[d]->apply(iv, d, map, value, this, &nr_m[0]);
        }
    }
    return applied;
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::applyBoundary(const AmrIntVect_t& iv,
                                                              const lo_t& dir,
                                                              umap_t& map,
                                                              const scalar_t& value)
{
    bc_mp[dir]->apply(iv, dir, map, value, this, &nr_m[0]);
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildLevelMask_m() {
    amrex::Periodicity period(AmrIntVect_t(D_DECL(0, 0, 0)));
    mask.reset(new mask_t(grids, dmap, 1, 1));
    mask->BuildMask(geom.Domain(), period,
                    Mask::COVERED, Mask::BNDRY,
                    Mask::PHYSBNDRY, Mask::INTERIOR);
    mask->FillBoundary(period);
}


template <class MatrixType, class VectorType>
const amr::AmrIntVect_t& AmrMultiGridLevel<MatrixType, VectorType>::refinement() const {
    return rr_m;
}


template <class MatrixType, class VectorType>
const amr::scalar_t* AmrMultiGridLevel<MatrixType, VectorType>::cellSize() const {
    return dx_m;
}


template <class MatrixType, class VectorType>
const amr::scalar_t& AmrMultiGridLevel<MatrixType, VectorType>::cellSize(lo_t dir) const {
    return dx_m[dir];
}


template <class MatrixType, class VectorType>
const amr::scalar_t* AmrMultiGridLevel<MatrixType, VectorType>::invCellSize() const {
    return invdx_m;
}


template <class MatrixType, class VectorType>
const amr::scalar_t& AmrMultiGridLevel<MatrixType, VectorType>::invCellSize(lo_t dir) const {
    return invdx_m[dir];
}


template <class MatrixType, class VectorType>
bool AmrMultiGridLevel<MatrixType, VectorType>::isValid(const AmrIntVect_t& iv) const {
    return ( iv.allGT(AmrIntVect_t(D_DECL(-1, -1, -1))) &&
             iv.allLT(AmrIntVect_t(D_DECL(nr_m[0], nr_m[1], nr_m[2]))) );
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildMap_m(const Teuchos::RCP<comm_t>& comm,
                                                           const Teuchos::RCP<node_t>& node)
{
    
    go_t localNumElements = 0;
    
    Teuchos::Array<go_t> globalindices;
    
    for (amrex::MFIter mfi(grids, dmap, true); mfi.isValid(); ++mfi) {
        const amrex::Box&    tbx  = mfi.tilebox();
        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));

                    go_t globalidx = serialize(iv);
                    
                    globalindices.push_back(globalidx);
                    
                    ++localNumElements;
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    /*
     * create map that specifies which processor gets which data
     */
    
    // get smallest global index of this level
    amrex::Box bx = grids.minimalBox();
    const int* lo = bx.loVect();
    AmrIntVect_t lowcorner(D_DECL(lo[0], lo[1], lo[2]));
    
    // where to start indexing
    go_t baseIndex = serialize(lowcorner);
    
    // numGlobalElements == N
    go_t N = grids.numPts();
    
    map_p = Teuchos::rcp( new dmap_t(N, globalindices, baseIndex, comm, node) );
}
