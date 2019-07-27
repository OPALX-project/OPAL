#include "AmrMultiGrid.h"

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include <AMReX_ParallelDescriptor.H>

#include <boost/filesystem.hpp>

#if AMR_MG_WRITE
    #include <iomanip>
#endif

AmrMultiGrid::AmrMultiGrid(AmrBoxLib* itsAmrObject_p,
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
                           const std::string& norm)
    : AmrPoissonSolver<AmrBoxLib>(itsAmrObject_p)
    , comm_mp( new comm_t( amrex::ParallelDescriptor::Communicator() ) )
    , nIter_m(0)
    , bIter_m(0)
    , maxiter_m(100)
    , nSweeps_m(nSweeps)
    , mglevel_m(0)
    , lbase_m(0)
    , lfine_m(0)
    , nlevel_m(1)
    , nBcPoints_m(0)
    , eps_m(1.0e-10)
    , verbose_m(false)
    , fname_m(OpalData::getInstance()->getInputBasename() + std::string(".solver"))
    , flag_m(std::ios::out)
{
    node_mp = KokkosClassic::Details::getNode<amr::node_t>(); //KokkosClassic::DefaultNode::getDefaultNode();

#if AMR_MG_TIMER
    this->initTimer_m();
#endif

    const Boundary bcs[AMREX_SPACEDIM] = {
        this->convertToEnumBoundary_m(bcx),
        this->convertToEnumBoundary_m(bcy),
        this->convertToEnumBoundary_m(bcz)
    };

    this->initPhysicalBoundary_m(&bcs[0]);

    smootherType_m = this->convertToEnumSmoother_m(smoother);

    norm_m = this->convertToEnumNorm_m(norm);

    const Interpolater interpolater = this->convertToEnumInterpolater_m(interp);
    this->initInterpolater_m(interpolater);

    // interpolater for crse-fine-interface
    this->initCrseFineInterp_m(Interpolater::LAGRANGE);

    // preconditioner
    const Preconditioner precond = this->convertToEnumPreconditioner_m(prec);
    this->initPrec_m(precond, rebalance, reuse);

    // base level solver
    const BaseSolver solver = this->convertToEnumBaseSolver_m(bsolver);
    this->initBaseSolver_m(solver, rebalance, reuse);

    if (boost::filesystem::exists(fname_m)) {
        flag_m = std::ios::app;
        INFOMSG("Appending solver information to existing file: " << fname_m << endl);
    } else {
        INFOMSG("Creating new file for solver information: " << fname_m << endl);
    }
}


void AmrMultiGrid::solve(AmrScalarFieldContainer_t &rho,
                         AmrScalarFieldContainer_t &phi,
                         AmrVectorFieldContainer_t &efield,
                         unsigned short baseLevel,
                         unsigned short finestLevel,
                         bool prevAsGuess)
{
    lbase_m = baseLevel;
    lfine_m = finestLevel;
    nlevel_m = lfine_m - lbase_m + 1;

    /* we cannot use the previous solution
     * if we have to regrid (AmrPoissonSolver::hasToRegrid())
     *
     * regrid_m is set in AmrBoxlib::regrid()
     */
    bool reset = !prevAsGuess;

    if ( this->regrid_m )
        reset = true;

    this->initLevels_m(rho, itsAmrObject_mp->Geom(), this->regrid_m);

    // build all necessary matrices and vectors
    this->setup_m(rho, phi, this->regrid_m);

    this->initGuess_m(reset);

    // actual solve
    scalar_t error = this->iterate_m();

    for (int lev = nlevel_m - 1; lev > -1; --lev) {
        averageDown_m(lev);
    }

    // write efield to AMReX
    this->computeEfield_m(efield);
    
    // copy solution back
    for (int lev = 0; lev < nlevel_m; ++lev) {
        int ilev = lbase_m + lev;

        phi[ilev]->setVal(0.0, phi[ilev]->nGrow());

        this->trilinos2amrex_m(lev, 0, *phi[ilev], mglevel_m[lev]->phi_p);
    }

    if ( verbose_m )
        this->writeSDDSData_m(error);

    // we can now reset
    this->regrid_m = false;
}


void AmrMultiGrid::setNumberOfSweeps(const std::size_t& nSweeps) {
    if ( nSweeps < 0 )
        throw OpalException("AmrMultiGrid::setNumberOfSweeps()",
                            "The number of smoothing sweeps needs to be non-negative!");

    nSweeps_m = nSweeps;
}


void AmrMultiGrid::setMaxNumberOfIterations(const std::size_t& maxiter) {
    if ( maxiter < 1 )
        throw OpalException("AmrMultiGrid::setMaxNumberOfIterations()",
                            "The max. number of iterations needs to be positive!");

    maxiter_m = maxiter;
}


std::size_t AmrMultiGrid::getNumIters() {
    return nIter_m;
}


AmrMultiGrid::scalar_t AmrMultiGrid::getLevelResidualNorm(lo_t level) {
    return evalNorm_m(mglevel_m[level]->residual_p);
}


void AmrMultiGrid::setVerbose(bool verbose) {
    verbose_m = verbose;
}


void AmrMultiGrid::initPhysicalBoundary_m(const Boundary* bc)
{
    // make sure it's reset
    nBcPoints_m = 0;

    for (unsigned int i = 0; i < AMREX_SPACEDIM; ++i) {
        switch ( bc[i] ) {
            case Boundary::DIRICHLET:
                bc_m[i].reset( new AmrDirichletBoundary<AmrMultiGridLevel_t>() );
                break;
            case Boundary::OPEN:
                bc_m[i].reset( new AmrOpenBoundary<AmrMultiGridLevel_t>() );
                break;
            case Boundary::PERIODIC:
                bc_m[i].reset( new AmrPeriodicBoundary<AmrMultiGridLevel_t>() );
                break;
            default:
                throw OpalException("AmrMultiGrid::initPhysicalBoundary_m()",
                                    "This type of boundary is not supported");
        }
        // we use the maximum in order to build matrices
        go_t tmp = bc_m[i]->getNumberOfPoints();
        if ( nBcPoints_m < tmp )
            nBcPoints_m = tmp;
    }
}


void AmrMultiGrid::initLevels_m(const amrex::Vector<AmrField_u>& rho,
                                const amrex::Vector<AmrGeometry_t>& geom,
                                bool regrid)
{
    if ( !regrid )
        return;

    // although we do a resize afterwards, we do this to be safe
    for (int lev = nlevel_m; lev < (int)mglevel_m.size(); ++lev) {
        mglevel_m[lev].reset(nullptr);
    }

    mglevel_m.resize(nlevel_m);

    amrex::Periodicity period(AmrIntVect_t(D_DECL(0, 0, 0)));

    AmrIntVect_t rr = AmrIntVect_t(D_DECL(2, 2, 2));

    for (int lev = 0; lev < nlevel_m; ++lev) {
        int ilev = lbase_m + lev;

        // do not initialize base level every time
        if (mglevel_m[lev] == nullptr || lev > lbase_m) {
            mglevel_m[lev].reset(new AmrMultiGridLevel_t(itsAmrObject_mp->getMeshScaling(),
                                                         rho[ilev]->boxArray(),
                                                         rho[ilev]->DistributionMap(),
                                                         geom[ilev],
                                                         rr,
                                                         bc_m,
                                                         comm_mp,
                                                         node_mp));
        } else {
            mglevel_m[lev]->buildLevelMask();
        }

        mglevel_m[lev]->refmask.reset(
            new AmrMultiGridLevel_t::mask_t(mglevel_m[lev]->grids,
                                            mglevel_m[lev]->dmap, 1, 2)
            );
        mglevel_m[lev]->refmask->setVal(AmrMultiGridLevel_t::Refined::NO, 2);
        mglevel_m[lev]->refmask->FillBoundary(period);

        amrex::BoxArray ba = mglevel_m[lev]->grids;
        ba.coarsen(rr);
        mglevel_m[lev]->crsemask.reset(
            new AmrMultiGridLevel_t::mask_t(ba,
                                            mglevel_m[lev]->dmap, 1, 2)
            );
        mglevel_m[lev]->crsemask->setVal(AmrMultiGridLevel_t::Refined::NO, 2);
    }

    for (int lev = 1; lev < nlevel_m; ++lev) {
        mglevel_m[lev]->crsemask->setVal(AmrMultiGridLevel_t::Refined::YES, 0);

        // used for boundary interpolation --> replaces expensive calls to isBoundary
        mglevel_m[lev]->crsemask->setDomainBndry(AmrMultiGridLevel_t::Mask::PHYSBNDRY,
                                                 mglevel_m[lev-1]->geom); //FIXME: geometry of lev - 1
        // really needed ?
        mglevel_m[lev]->crsemask->FillBoundary(period);
    }

    /* to complete initialization we need to fill
     * the mask of refinement
     */
    for (int lev = 0; lev < nlevel_m-1; ++lev) {
        // get boxarray with refined cells
        amrex::BoxArray ba = mglevel_m[lev]->grids;
        ba.refine(rr);
        ba = amrex::intersect(mglevel_m[lev+1]->grids, ba);
        ba.coarsen(rr);

        // refined cells
        amrex::DistributionMapping dmap(ba, comm_mp->getSize());
        AmrMultiGridLevel_t::mask_t refined(ba, dmap, 1, 0);
        refined.setVal(AmrMultiGridLevel_t::Refined::YES);
//      refined.setDomainBndry(AmrMultiGridLevel_t::Mask::PHYSBNDRY, mglevel_m[lev]->geom);

        // fill intersection with YES
        mglevel_m[lev]->refmask->copy(refined, 0, 0, 1, 0, 2);

        /* physical boundary cells will never be refined cells
         * since they are ghost cells
         */
        mglevel_m[lev]->refmask->setDomainBndry(AmrMultiGridLevel_t::Mask::PHYSBNDRY,
                                                mglevel_m[lev]->geom);

        mglevel_m[lev]->refmask->FillBoundary(period);
    }
}


void AmrMultiGrid::clearMasks_m() {
    for (int lev = 0; lev < nlevel_m; ++lev) {
        mglevel_m[lev]->refmask.reset(nullptr);
        mglevel_m[lev]->crsemask.reset(nullptr);
        mglevel_m[lev]->mask.reset(nullptr);
    }
}


void AmrMultiGrid::initGuess_m(bool reset) {
    if ( !reset )
        return;

    // reset
    for (int lev = 0; lev < nlevel_m; ++lev)
        mglevel_m[lev]->phi_p->putScalar(0.0);
}


AmrMultiGrid::scalar_t AmrMultiGrid::iterate_m() {

    // initial error
    std::vector<scalar_t> rhsNorms;
    std::vector<scalar_t> resNorms;

    this->initResidual_m(rhsNorms, resNorms);

    std::for_each(rhsNorms.begin(), rhsNorms.end(),
                  [this](double& val){ val *= eps_m; });

    nIter_m = 0;
    bIter_m = 0;

    while ( !isConverged_m(rhsNorms, resNorms) && nIter_m < maxiter_m ) {

        relax_m(lfine_m);

//         /* in contrast to algorithm, we average down now
//          * --> potential is valid also on coarse covered
//          * cells
//          * --> however, it may take 1-2 iterations longer
//          */
//         for (int lev = nlevel_m - 1; lev > -1; --lev) {
//             averageDown_m(lev);
//         }

        // update residual
        for (int lev = 0; lev < nlevel_m; ++lev) {

            this->residual_m(lev,
                             mglevel_m[lev]->residual_p,
                             mglevel_m[lev]->rho_p,
                             mglevel_m[lev]->phi_p);
        }


        for (lo_t lev = 0; lev < nlevel_m; ++lev)
            resNorms[lev] = getLevelResidualNorm(lev);

        ++nIter_m;

#if AMR_MG_WRITE
        this->writeResidualNorm_m();
#endif

        bIter_m += solver_mp->getNumIters();
    }

    return std::accumulate(resNorms.begin(),
                           resNorms.end(), 0.0,
                           std::plus<scalar_t>());
}


bool AmrMultiGrid::isConverged_m(std::vector<scalar_t>& rhsNorms,
                                 std::vector<scalar_t>& resNorms)
{
    return std::equal(resNorms.begin(), resNorms.end(),
                      rhsNorms.begin(),
                      std::less<scalar_t>());
}


void AmrMultiGrid::residual_m(const lo_t& level,
                              Teuchos::RCP<vector_t>& r,
                              const Teuchos::RCP<vector_t>& b,
                              const Teuchos::RCP<vector_t>& x)
{
    /*
     * r = b - A*x
     */
    if ( level < lfine_m ) {

        vector_t fine2crse(mglevel_m[level]->Awf_p->getDomainMap(), true);

        // get boundary for
        if ( mglevel_m[level]->Bfine_p != Teuchos::null ) {
            mglevel_m[level]->Bfine_p->apply(*mglevel_m[level+1]->phi_p, fine2crse);
        }

        // operation: fine2crse += A * x
        mglevel_m[level]->Awf_p->apply(*x, fine2crse,
                                       Teuchos::NO_TRANS,
                                       scalar_t(1.0),
                                       scalar_t(1.0));

        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            // operation: fine2crse += B * phi^(l-1)
            mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p,
                                             fine2crse,
                                             Teuchos::NO_TRANS,
                                             scalar_t(1.0),
                                             scalar_t(1.0));
        }

        Teuchos::RCP<vector_t> uncov_Ax = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
        mglevel_m[level]->UnCovered_p->apply(fine2crse, *uncov_Ax);

        Teuchos::RCP<vector_t> uncov_b = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );

        mglevel_m[level]->UnCovered_p->apply(*b, *uncov_b);

        // ONLY subtract coarse rho
//         mglevel_m[level]->residual_p->putScalar(0.0);

        r->update(1.0, *uncov_b, -1.0, *uncov_Ax, 0.0);

    } else {
        /* finest level: Awf_p == Anf_p
         */
        Teuchos::RCP<vector_t> Ax = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
        mglevel_m[level]->Anf_p->apply(*x, *Ax);
        
        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            // operationr: Ax += B * phi^(l-1)
            mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p,
                                             *Ax,
                                             Teuchos::NO_TRANS,
                                             scalar_t(1.0),
                                             scalar_t(1.0));
        }
        r->update(1.0, *b, -1.0, *Ax, 0.0);
    }
}


void AmrMultiGrid::relax_m(const lo_t& level) {

    if ( level == lfine_m ) {

        if ( level == lbase_m ) {
            /* Anf_p == Awf_p
             */
            Teuchos::RCP<vector_t> Ax = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
            mglevel_m[level]->Anf_p->apply(*mglevel_m[level]->phi_p, *Ax);
            mglevel_m[level]->residual_p->update(1.0, *mglevel_m[level]->rho_p, -1.0, *Ax, 0.0);
            
        } else {
            this->residual_no_fine_m(level,
                                     mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     mglevel_m[level-1]->phi_p,
                                     mglevel_m[level]->rho_p);
        }
    }

    if ( level > 0 ) {
        // phi^(l, save) = phi^(l)
        Teuchos::RCP<vector_t> phi_save = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        Tpetra::deep_copy(*phi_save, *mglevel_m[level]->phi_p);

        mglevel_m[level-1]->error_p->putScalar(0.0);

        // smoothing
        this->smooth_m(level,
                       mglevel_m[level]->error_p,
                       mglevel_m[level]->residual_p);


        // phi = phi + e
        mglevel_m[level]->phi_p->update(1.0, *mglevel_m[level]->error_p, 1.0);

        /*
         * restrict
         */
        this->restrict_m(level);

        this->relax_m(level - 1);

        /*
         * prolongate / interpolate
         */
        this->prolongate_m(level);

        // residual update
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        this->residual_no_fine_m(level, tmp,
                                 mglevel_m[level]->error_p,
                                 mglevel_m[level-1]->error_p,
                                 mglevel_m[level]->residual_p);

        Tpetra::deep_copy(*mglevel_m[level]->residual_p, *tmp);

        // delta error
        Teuchos::RCP<vector_t> derror = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );

        // smoothing
        this->smooth_m(level, derror, mglevel_m[level]->residual_p);

        // e^(l) += de^(l)
        mglevel_m[level]->error_p->update(1.0, *derror, 1.0);

        // phi^(l) = phi^(l, save) + e^(l)
        mglevel_m[level]->phi_p->update(1.0, *phi_save, 1.0, *mglevel_m[level]->error_p, 0.0);

    } else {
        // e = A^(-1)r
#if AMR_MG_TIMER
    IpplTimings::startTimer(bottomTimer_m);
#endif

        solver_mp->solve(mglevel_m[level]->error_p,
                         mglevel_m[level]->residual_p);

#if AMR_MG_TIMER
    IpplTimings::stopTimer(bottomTimer_m);
#endif
        // phi = phi + e
        mglevel_m[level]->phi_p->update(1.0, *mglevel_m[level]->error_p, 1.0);
    }
}


void AmrMultiGrid::residual_no_fine_m(const lo_t& level,
                                      Teuchos::RCP<vector_t>& result,
                                      const Teuchos::RCP<vector_t>& rhs,
                                      const Teuchos::RCP<vector_t>& crs_rhs,
                                      const Teuchos::RCP<vector_t>& b)
{
    vector_t crse2fine(mglevel_m[level]->Anf_p->getDomainMap(), true);

    // get boundary for
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->apply(*crs_rhs, crse2fine);
    }

    // operation: crse2fine = 1.0 * crse2fine  + 1.0 * A^(l) * rhs
    mglevel_m[level]->Anf_p->apply(*rhs, crse2fine,
                                   Teuchos::NO_TRANS,
                                   scalar_t(1.0),
                                   scalar_t(1.0));

    result->update(1.0, *b, -1.0, crse2fine, 0.0);
}


#if AMR_MG_WRITE
void AmrMultiGrid::writeResidualNorm_m() {
    scalar_t err = 0.0;

    std::ofstream out;
    if ( comm_mp->getRank() == 0 )
        out.open("residual.dat", std::ios::app);

    for (int lev = 0; lev < nlevel_m; ++lev) {
        scalar_t tmp = evalNorm_m(mglevel_m[lev]->residual_p);

        if ( comm_mp->getRank() == 0 )
            out << std::setw(15) << std::right << tmp;
    }

    if ( comm_mp->getRank() == 0 )
        out.close();
}
#endif


typename AmrMultiGrid::scalar_t
AmrMultiGrid::evalNorm_m(const Teuchos::RCP<const vector_t>& x)
{
    scalar_t norm = 0.0;

    switch ( norm_m ) {
        case Norm::L1:
        {
            norm = x->norm1();
            break;
        }
        case Norm::L2:
        {
            norm = x->norm2();
            break;
        }
        case Norm::LINF:
        {
            norm = x->normInf();
            break;
        }
        default:
            throw OpalException("AmrMultiGrid::evalNorm_m()",
                                "This type of norm not suppported.");
    }
    return norm;
}


void AmrMultiGrid::initResidual_m(std::vector<scalar_t>& rhsNorms,
                                  std::vector<scalar_t>& resNorms)
{
    rhsNorms.clear();
    resNorms.clear();

#if AMR_MG_WRITE
    std::ofstream out;

    if ( comm_mp->getRank() == 0) {
        out.open("residual.dat", std::ios::out);

        for (int lev = 0; lev < nlevel_m; ++lev)
            out << std::setw(14) << std::right << "level" << lev;
        out << std::endl;
    }
#endif

    for (int lev = 0; lev < nlevel_m; ++lev) {
        this->residual_m(lev,
                         mglevel_m[lev]->residual_p,
                         mglevel_m[lev]->rho_p,
                         mglevel_m[lev]->phi_p);

        resNorms.push_back(evalNorm_m(mglevel_m[lev]->residual_p));

#if AMR_MG_WRITE
        if ( comm_mp->getRank() == 0 )
            out << std::setw(15) << std::right << resNorms.back();
#endif

        rhsNorms.push_back(evalNorm_m(mglevel_m[lev]->rho_p));
    }

#if AMR_MG_WRITE
    if ( comm_mp->getRank() == 0 )
        out.close();
#endif
}


void AmrMultiGrid::computeEfield_m(AmrVectorFieldContainer_t& efield) {
#if AMR_MG_TIMER
    IpplTimings::startTimer(efieldTimer_m);
#endif
    Teuchos::RCP<vector_t> efield_p = Teuchos::null;
    for (int lev = nlevel_m - 1; lev > -1; --lev) {
        int ilev = lbase_m + lev;

        efield_p = Teuchos::rcp( new vector_t(mglevel_m[lev]->map_p, false) );


        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            efield[ilev][d]->setVal(0.0, efield[ilev][d]->nGrow());
            efield_p->putScalar(0.0);
            mglevel_m[lev]->G_p[d]->apply(*mglevel_m[lev]->phi_p, *efield_p);
            this->trilinos2amrex_m(lev, 0, *efield[ilev][d], efield_p);
        }
    }
#if AMR_MG_TIMER
    IpplTimings::stopTimer(efieldTimer_m);
#endif
}


void AmrMultiGrid::setup_m(const amrex::Vector<AmrField_u>& rho,
                           const amrex::Vector<AmrField_u>& phi,
                           const bool& matrices)
{
#if AMR_MG_TIMER
    IpplTimings::startTimer(buildTimer_m);
#endif

    if ( lbase_m == lfine_m )
        this->buildSingleLevel_m(rho, phi, matrices);
    else
        this->buildMultiLevel_m(rho, phi, matrices);

    mglevel_m[lfine_m]->error_p->putScalar(0.0);

    if ( matrices ) {
        this->clearMasks_m();
        // set the bottom solve operator
        if (!solver_mp->hasOperator()) {
            solver_mp->setOperator(mglevel_m[lbase_m]->Anf_p, mglevel_m[0].get());
        }
    }

#if AMR_MG_TIMER
    IpplTimings::stopTimer(buildTimer_m);
#endif
}


void AmrMultiGrid::buildSingleLevel_m(const amrex::Vector<AmrField_u>& rho,
                                      const amrex::Vector<AmrField_u>& phi,
                                      const bool& matrices)
{
    this->open_m(lbase_m, matrices);

    const scalar_t* invdx = mglevel_m[lbase_m]->invCellSize();

    const scalar_t invdx2[] = {
        D_DECL( invdx[0] * invdx[0],
                invdx[1] * invdx[1],
                invdx[2] * invdx[2] )
    };

    if (matrices) {
        for (amrex::MFIter mfi(*mglevel_m[lbase_m]->mask, true);
             mfi.isValid(); ++mfi)
        {
            const box_t&       tbx = mfi.tilebox();
            const basefab_t&  mfab = (*mglevel_m[lbase_m]->mask)[mfi];
            const farraybox_t& rhofab = (*rho[lbase_m])[mfi];
            const farraybox_t& pfab = (*phi[lbase_m])[mfi];

            const int* lo = tbx.loVect();
            const int* hi = tbx.hiVect();

            for (int i = lo[0]; i <= hi[0]; ++i) {
                for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                    for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                        AmrIntVect_t iv(D_DECL(i, j, k));
                        go_t gidx = mglevel_m[lbase_m]->serialize(iv);

                        if (!solver_mp->hasOperator()) {
                            this->buildNoFinePoissonMatrix_m(lbase_m, gidx, iv, mfab, invdx2);
                            this->buildGradientMatrix_m(lbase_m, gidx, iv, mfab, invdx);
                        }

                        mglevel_m[lbase_m]->rho_p->replaceGlobalValue(gidx, rhofab(iv, 0));
                        mglevel_m[lbase_m]->phi_p->replaceGlobalValue(gidx, pfab(iv, 0));

#if AMREX_SPACEDIM == 3
                    }
#endif
                }
            }
        }
    } else {
        this->buildDensityVector_m(lbase_m, *rho[lbase_m]);
    }

    this->close_m(lbase_m, matrices);

    if (matrices) {
        mglevel_m[lbase_m]->Awf_p = Teuchos::null;
        mglevel_m[lbase_m]->UnCovered_p = Teuchos::null;
    }
}


void AmrMultiGrid::buildMultiLevel_m(const amrex::Vector<AmrField_u>& rho,
                                     const amrex::Vector<AmrField_u>& phi,
                                     const bool& matrices)
{
    // the base level has no smoother --> nlevel_m - 1
    if ( matrices ) {
        // although we do a resize afterwards, we do this to be safe
        for (int lev = nlevel_m-1; lev < (int)smoother_m.size(); ++lev) {
            smoother_m[lev].reset();
        }
        smoother_m.resize(nlevel_m-1);
    }
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        this->open_m(lev, matrices);

        int ilev = lbase_m + lev;

        // find all coarse cells that are covered by fine cells
//         AmrIntVect_t rr = mglevel_m[lev]->refinement();

        const scalar_t* invdx = mglevel_m[lev]->invCellSize();

        const scalar_t invdx2[] = {
            D_DECL( invdx[0] * invdx[0],
                    invdx[1] * invdx[1],
                    invdx[2] * invdx[2] )
        };

        if ( matrices ) {
            for (amrex::MFIter mfi(*mglevel_m[lev]->mask, true);
                 mfi.isValid(); ++mfi)
            {
                const box_t&       tbx = mfi.tilebox();
                const basefab_t&  mfab = (*mglevel_m[lev]->mask)[mfi];
                const basefab_t&  rfab = (*mglevel_m[lev]->refmask)[mfi];
                const basefab_t&  cfab = (*mglevel_m[lev]->crsemask)[mfi];
                const farraybox_t& rhofab = (*rho[ilev])[mfi];
                const farraybox_t& pfab = (*phi[ilev])[mfi];

                const int* lo = tbx.loVect();
                const int* hi = tbx.hiVect();

                for (int i = lo[0]; i <= hi[0]; ++i) {
                    int ii = i << 1;
                    for (int j = lo[1]; j <= hi[1]; ++j) {
                        int jj = j << 1;
#if AMREX_SPACEDIM == 3
                        for (int k = lo[2]; k <= hi[2]; ++k) {
                            int kk = k << 1;
#endif
                            AmrIntVect_t iv(D_DECL(i, j, k));
                            go_t gidx = mglevel_m[lev]->serialize(iv);

                            this->buildRestrictionMatrix_m(lev, gidx, iv,
                                                           D_DECL(ii, jj, kk), rfab);

                            this->buildInterpolationMatrix_m(lev, gidx, iv, cfab);

                            this->buildCrseBoundaryMatrix_m(lev, gidx, iv, mfab,
                                                            cfab, invdx2);

                            this->buildFineBoundaryMatrix_m(lev, gidx, iv,
                                                            mfab, rfab, cfab);

                            this->buildCompositePoissonMatrix_m(lev, gidx, iv, mfab,
                                                                rfab, cfab, invdx2);

                            if (lev > lbase_m || (lev == lbase_m && !solver_mp->hasOperator())) {
                                this->buildNoFinePoissonMatrix_m(lev, gidx, iv, mfab, invdx2);
                                this->buildGradientMatrix_m(lev, gidx, iv, mfab, invdx);
                            }

                            mglevel_m[lev]->rho_p->replaceGlobalValue(gidx, rhofab(iv, 0));
                            mglevel_m[lev]->phi_p->replaceGlobalValue(gidx, pfab(iv, 0));
#if AMREX_SPACEDIM == 3
                        }
#endif
                    }
                }
            }
        } else {
            for (lo_t lev = 0; lev < nlevel_m; ++lev) {
                int ilev = lbase_m + lev;
                this->buildDensityVector_m(lev, *rho[ilev]);
            }
        }

        this->close_m(lev, matrices);

        if ( matrices && lev > lbase_m ) {
            smoother_m[lev-1].reset( new AmrSmoother(mglevel_m[lev]->Anf_p,
                                                     smootherType_m, nSweeps_m) );
        }
    }
}


void AmrMultiGrid::open_m(const lo_t& level,
                          const bool& matrices)
{
    if ( matrices ) {

        if ( level > lbase_m ) {

            /*
             * interpolation matrix
             */

            int nNeighbours = (nBcPoints_m + 1) * interp_mp->getNumberOfPoints();

            mglevel_m[level]->I_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                               nNeighbours,
                                                               Tpetra::StaticProfile) );

            /*
             * coarse boundary matrix
             */

            nNeighbours = 2 * AMREX_SPACEDIM * nBcPoints_m *
                          2 * AMREX_SPACEDIM * interface_mp->getNumberOfPoints();

            mglevel_m[level]->Bcrse_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );

        }


        if ( level < lfine_m ) {

            /*
             * restriction matrix
             */

            // refinement 2
            int nNeighbours = AMREX_D_TERM(2, * 2, * 2);

            mglevel_m[level]->R_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );

            /*
             * fine boundary matrix
             */

            //                              refinement 2
            nNeighbours = 2 * AMREX_SPACEDIM * AMREX_D_TERM(2, * 2, * 2);

            mglevel_m[level]->Bfine_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );

        }

        /*
         * no-fine Poisson matrix
         */

        int nPhysBoundary = 2 * AMREX_SPACEDIM * nBcPoints_m;

        // number of internal stencil points
        int nIntBoundary = AMREX_SPACEDIM * interface_mp->getNumberOfPoints();

        int nEntries = (AMREX_SPACEDIM << 1) + 2 /* plus boundaries */ + nPhysBoundary + nIntBoundary;

        if (level > lbase_m || (level == lbase_m && !solver_mp->hasOperator())) {
            mglevel_m[level]->Anf_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nEntries,
                             Tpetra::StaticProfile) );
        }

        /*
         * with-fine / composite Poisson matrix
         */
        if ( lbase_m != lfine_m ) {
            nEntries = (AMREX_SPACEDIM << 1) + 5 /* plus boundaries */ + nPhysBoundary + nIntBoundary;
    
            mglevel_m[level]->Awf_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nEntries,
                             Tpetra::StaticProfile) );
        
            /*
             * uncovered cells matrix
             */
            mglevel_m[level]->UnCovered_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, 1,
                             Tpetra::StaticProfile) );
        }
        
        /*
         * gradient matrices
         */
        nEntries = 11;

        if (level > lbase_m || (level == lbase_m && !solver_mp->hasOperator())) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                mglevel_m[level]->G_p[d] = Teuchos::rcp(
                    new matrix_t(mglevel_m[level]->map_p, nEntries,
                                 Tpetra::StaticProfile) );
            }
        }
    }

    mglevel_m[level]->rho_p = Teuchos::rcp(
        new vector_t(mglevel_m[level]->map_p, false) );

    if ( matrices ) {
        mglevel_m[level]->phi_p = Teuchos::rcp(
            new vector_t(mglevel_m[level]->map_p, false) );
    }
}


void AmrMultiGrid::close_m(const lo_t& level,
                           const bool& matrices)
{
    if ( matrices ) {
        if ( level > lbase_m ) {

            mglevel_m[level]->I_p->fillComplete(mglevel_m[level-1]->map_p,  // col map (domain map)
                                                mglevel_m[level]->map_p);   // row map (range map)

            mglevel_m[level]->Bcrse_p->fillComplete(mglevel_m[level-1]->map_p,  // col map
                                                    mglevel_m[level]->map_p);   // row map
        }

        if ( level < lfine_m ) {

            mglevel_m[level]->R_p->fillComplete(mglevel_m[level+1]->map_p,
                                                  mglevel_m[level]->map_p);

            mglevel_m[level]->Bfine_p->fillComplete(mglevel_m[level+1]->map_p,
                                                    mglevel_m[level]->map_p);
        }

        if (level > lbase_m || (level == lbase_m && !solver_mp->hasOperator())) {
            mglevel_m[level]->Anf_p->fillComplete();

            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                mglevel_m[level]->G_p[d]->fillComplete();
        }

        if ( lbase_m != lfine_m ) {
            mglevel_m[level]->Awf_p->fillComplete();

            mglevel_m[level]->UnCovered_p->fillComplete();
        }
    }
}


void AmrMultiGrid::buildNoFinePoissonMatrix_m(const lo_t& level,
                                              const go_t& gidx,
                                              const AmrIntVect_t& iv,
                                              const basefab_t& mfab,
                                              const scalar_t* invdx2)
{
    /*
     * Laplacian of "no fine"
     */

    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */

    umap_t map;
    indices_t indices;
    coefficients_t values;

    /*
     * check neighbours in all directions (Laplacian stencil --> cross)
     */
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            AmrIntVect_t biv = iv;
            biv[d] += shift;

            switch ( mfab(biv) )
            {
                case AmrMultiGridLevel_t::Mask::INTERIOR:
                case AmrMultiGridLevel_t::Mask::COVERED: // covered --> interior cell
                {
                    map[mglevel_m[level]->serialize(biv)] += invdx2[d];
                    break;
                }
                case AmrMultiGridLevel_t::Mask::BNDRY:
                {
                    // boundary cell
                    // only level > 0 have this kind of boundary
#ifndef NDEBUG
                    if ( level == lbase_m )
                        throw OpalException("AmrMultiGrid::buildNoFinePoissonMatrix_m()",
                                            "Error in mask for level "
                                            + std::to_string(level) + "!");
#endif
                    /* Dirichlet boundary conditions from coarser level.
                     */
                    interface_mp->fine(biv, map, invdx2[d], d, -shift,
                                       mglevel_m[level].get());
                    break;
                }
                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                {
                    // physical boundary cell
                    mglevel_m[level]->applyBoundary(biv, d, map,
                                                    invdx2[d] /*matrix coefficient*/);
                    break;
                }
                default:
                    throw OpalException("AmrMultiGrid::buildNoFinePoissonMatrix_m()",
                                        "Error in mask for level "
                                        + std::to_string(level) + "!");
            }
        }
    }

    // check center
    map[gidx] += AMREX_D_TERM(- 2.0 * invdx2[0],
                              - 2.0 * invdx2[1],
                              - 2.0 * invdx2[2]);

    this->map2vector_m(map, indices, values);

    mglevel_m[level]->Anf_p->insertGlobalValues(gidx,
                                                indices.size(),
                                                &values[0],
                                                &indices[0]);
}


void AmrMultiGrid::buildCompositePoissonMatrix_m(const lo_t& level,
                                                 const go_t& gidx,
                                                 const AmrIntVect_t& iv,
                                                 const basefab_t& mfab,
                                                 const basefab_t& rfab,
                                                 const basefab_t& cfab,
                                                 const scalar_t* invdx2)
{
    /*
     * Laplacian of "with fine"
     *
     * For the finest level: Awf == Anf
     */
    if ( rfab(iv) == AmrMultiGridLevel_t::Refined::YES ) //|| lbase_m != lfine_m )
        return;
    /*
     * Only cells that are not refined
     */

    /*
     * 1D not supported by AmrMultiGrid
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */

    umap_t map;
    indices_t indices;
    coefficients_t values;

    /*
     * Only cells that are not refined
     */

    /*
     * check neighbours in all directions (Laplacian stencil --> cross)
     */
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            AmrIntVect_t biv = iv;
            biv[d] += shift;

            if ( rfab(biv) != AmrMultiGridLevel_t::Refined::YES )
            {
                /*
                 * It can't be a refined cell!
                 */
                switch ( mfab(biv) )
                {
                case AmrMultiGridLevel_t::Mask::COVERED:
                    // covered --> interior cell
                case AmrMultiGridLevel_t::Mask::INTERIOR:
                {
                    map[mglevel_m[level]->serialize(biv)] += invdx2[d];
                    map[gidx] -= invdx2[d]; // add center once
                    break;
                }
                case AmrMultiGridLevel_t::Mask::BNDRY:
                {
                    // boundary cell
                    // only level > 0 have this kind of boundary
#ifndef NDEBUG
                    if ( level == lbase_m )
                        throw OpalException("AmrMultiGrid::buildCompositePoissonMatrix_m()",
                                            "Error in mask for level "
                                            + std::to_string(level) + "!");
#endif

                    /* We are on the fine side of the crse-fine interface
                     * --> normal stencil --> no flux matching required
                     * --> interpolation of fine ghost cell required
                     * (used together with Bcrse)
                     */

                    /* Dirichlet boundary conditions from coarser level.
                     */
                    interface_mp->fine(biv, map, invdx2[d], d, -shift,
                                       mglevel_m[level].get());

                    // add center once
                    map[gidx] -= invdx2[d];
                    break;
                }
                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                {
                    // physical boundary cell
                    mglevel_m[level]->applyBoundary(biv, d, map,
                                                    invdx2[d] /*matrix coefficient*/);

                    // add center once
                    map[gidx] -= invdx2[d];
                    break;
                }
                default:
                    throw OpalException("AmrMultiGrid::buildCompositePoissonMatrix_m()",
                                        "Error in mask for level "
                                        + std::to_string(level) + "!");
                }
            } else {
                /*
                 * If neighbour cell is refined, we are on the coarse
                 * side of the crse-fine interface --> flux matching
                 * required --> interpolation of fine ghost cell
                 * (used together with Bfine)
                 */


                // flux matching, coarse part

                /* 2D --> 2 fine cells to compute flux per coarse-fine-interace --> avg = 2
                 * 3D --> 4 fine cells to compute flux per coarse-fine-interace --> avg = 4
                 *
                 * @precondition: refinement of 2
                 */
                // top and bottom for all directions
                const scalar_t* invcdx = mglevel_m[level]->invCellSize();
                const scalar_t* invfdx = mglevel_m[level+1]->invCellSize();
                scalar_t invavg = AMREX_D_PICK(1.0, 0.5, 0.25);
                scalar_t value = - invavg * invcdx[d] * invfdx[d];

                for (int d1 = 0; d1 < 2; ++d1) {
#if AMREX_SPACEDIM == 3
                    for (int d2 = 0; d2 < 2; ++d2) {
#endif

                        /* in order to get a top iv --> needs to be odd value in "d"
                         * in order to get a bottom iv --> needs to be even value in "d"
                         */
                        AmrIntVect_t fake(D_DECL(0, 0, 0));

                        fake[(d+1)%AMREX_SPACEDIM] = d1;
#if AMREX_SPACEDIM == 3
                        fake[(d+2)%AMREX_SPACEDIM] = d2;
#endif
                        interface_mp->coarse(iv, map, value, d, shift, rfab,
                                             fake, mglevel_m[level].get());

#if AMREX_SPACEDIM == 3
                    }
#endif
                }
            }
        }
    }

    this->map2vector_m(map, indices, values);

    mglevel_m[level]->Awf_p->insertGlobalValues(gidx,
                                                indices.size(),
                                                &values[0],
                                                &indices[0]);

    scalar_t vv = 1.0;
    mglevel_m[level]->UnCovered_p->insertGlobalValues(gidx,
                                                      1,
                                                      &vv,
                                                      &gidx);
}


void AmrMultiGrid::buildRestrictionMatrix_m(const lo_t& level,
                                            const go_t& gidx,
                                            const AmrIntVect_t& iv,
                                            D_DECL(const go_t& ii,
                                                   const go_t& jj,
                                                   const go_t& kk),
                                            const basefab_t& rfab)
{
    /*
     * x^(l) = R * x^(l+1)
     */

    // finest level does not need to have a restriction matrix
    if ( rfab(iv) == AmrMultiGridLevel_t::Refined::NO || level == lfine_m )
        return;

    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     *
     *
     */
    indices_t indices;
    indices.reserve(2 << (AMREX_SPACEDIM - 1));
    coefficients_t values;
    values.reserve(2 << (AMREX_SPACEDIM -1));

    // neighbours
    for (int iref = 0; iref < 2; ++iref) {
        for (int jref = 0; jref < 2; ++jref) {
#if AMREX_SPACEDIM == 3
            for (int kref = 0; kref < 2; ++kref) {
#endif
                AmrIntVect_t riv(D_DECL(ii + iref, jj + jref, kk + kref));

                indices.push_back( mglevel_m[level+1]->serialize(riv) );
                values.push_back( AMREX_D_PICK(0.5, 0.25, 0.125) );
#if AMREX_SPACEDIM == 3
            }
#endif
        }
    }

    mglevel_m[level]->R_p->insertGlobalValues(gidx,
                                              indices.size(),
                                              &values[0],
                                              &indices[0]);
}


void AmrMultiGrid::buildInterpolationMatrix_m(const lo_t& level,
                                              const go_t& gidx,
                                              const AmrIntVect_t& iv,
                                              const basefab_t& cfab)
{
    /* crse: level - 1
     * fine (this): level
     */

    /*
     * This does not include ghost cells
     * --> no boundaries
     *
     * x^(l) = I * x^(l-1)
     */

    if ( level == lbase_m )
        return;

    umap_t map;
    indices_t indices;
    coefficients_t values;

    /*
     * we need boundary + indices from coarser level
     */
    interp_mp->stencil(iv, cfab,  map, 1.0, mglevel_m[level-1].get());

    this->map2vector_m(map, indices, values);

    mglevel_m[level]->I_p->insertGlobalValues(gidx,
                                              indices.size(),
                                              &values[0],
                                              &indices[0]);
}


void AmrMultiGrid::buildCrseBoundaryMatrix_m(const lo_t& level,
                                             const go_t& gidx,
                                             const AmrIntVect_t& iv,
                                             const basefab_t& mfab,
                                             const basefab_t& cfab,
                                             const scalar_t* invdx2)
{
    /*
     * fine (this): level
     * coarse:      level - 1
     */

    // the base level has only physical boundaries
    if ( level == lbase_m )
        return;

    // iv is a fine cell

    umap_t map;
    indices_t indices;
    coefficients_t values;

    // check its neighbours to see if at crse-fine interface
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            // neighbour
            AmrIntVect_t niv = iv;
            niv[d] += shift;

            if ( mfab(niv) == AmrMultiGridLevel_t::Mask::BNDRY )
            {
                // neighbour does not belong to fine grids
                // includes no cells at physical boundary

                // coarse cell that is not refined
                AmrIntVect_t civ = iv;
                civ[d] += shift;
                civ.coarsen(mglevel_m[level]->refinement());

                // we need boundary + indices from coarser level
                // we need normalization by mesh size squared (of fine cell size)
                // --> Laplacian for fine cell
                // (Dirichlet boundary for coarse --> fine)
                interface_mp->coarse(civ, map, invdx2[d], d, shift, cfab,
                                     iv, mglevel_m[level-1].get());
            }
#ifndef NDEBUG
            else if ( mfab(niv) == AmrMultiGridLevel_t::Mask::PHYSBNDRY ) {
                throw OpalException("AmrMultiGrid::buildCrseBoundaryMatrix_m()",
                                    "Fine meshes shouldn't be connected "
                                    "to physical (i.e. mesh) boundary!");
            }
#endif
        }
    }

    this->map2vector_m(map, indices, values);

    mglevel_m[level]->Bcrse_p->insertGlobalValues(gidx,
                                                  indices.size(),
                                                  &values[0],
                                                  &indices[0]);
}


void AmrMultiGrid::buildFineBoundaryMatrix_m(const lo_t& level,
                                             const go_t& gidx,
                                             const AmrIntVect_t& iv,
                                             const basefab_t& mfab,
                                             const basefab_t& rfab,
                                             const basefab_t& cfab)
{
    /* fine: level + 1
     * coarse (this): level
     */

    // the finest level does not need data from a finer level
    if ( rfab(iv) == AmrMultiGridLevel_t::Refined::YES || level == lfine_m )
        return;

    const scalar_t* invcdx = mglevel_m[level]->invCellSize();
    const scalar_t* invfdx = mglevel_m[level+1]->invCellSize();

    // inverse of number of fine cell gradients
    scalar_t invavg = AMREX_D_PICK(1, 0.5, 0.25);

    umap_t map;
    indices_t indices;
    coefficients_t values;

    auto fill = [&](umap_t& map,
                    D_DECL(int ii, int jj, int kk),
                    int* begin, int* end, int d,
                    const AmrIntVect_t& iv, int shift,
                    int sign)
    {
        for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
            for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
#if AMREX_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     */
                    AmrIntVect_t riv(D_DECL(iref, jref, kref));

                    if ( (riv[d] >> 1) /*refinement*/ == iv[d] ) {
                        /* the fine cell is on the coarse side --> fine
                         * ghost cell --> we need to interpolate
                         */
                        scalar_t value = - invavg * invcdx[d] * invfdx[d];

                        interface_mp->fine(riv, map, value, d, shift,
                                           mglevel_m[level+1].get());
                    } else {
                        scalar_t value = invavg * invcdx[d] * invfdx[d];
                        map[mglevel_m[level+1]->serialize(riv)] += value;
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    };


    /*
     * iv is a coarse cell that got not refined
     *
     * --> check all neighbours to see if at crse-fine
     * interface
     */
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            // neighbour
            AmrIntVect_t covered = iv;
            covered[d] += shift;

            if ( rfab(covered) == AmrMultiGridLevel_t::Refined::YES &&
                 mfab(covered) != AmrMultiGridLevel_t::PHYSBNDRY )
            {
                // neighbour is covered by fine cells

                /*
                 * "shift" is the amount to a coarse cell that got refined
                 * "d" is the direction to shift
                 *
                 * --> check all covered neighbour cells
                 */

                /* we need to iterate over correct fine cells. It depends
                 * on the orientation of the interface
                 */
                int begin[AMREX_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
                int end[AMREX_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };

                /*
                 * neighbour cell got refined but is not on physical boundary
                 * --> we are at a crse-fine interface
                 *
                 * we need now to find out which fine cells
                 * are required to satisfy the flux matching
                 * condition
                 */

                switch ( shift ) {
                    case -1:
                    {
                        // --> interface is on the lower face
                        int ii = iv[0] << 1; // refinemet in x
                        int jj = iv[1] << 1; // refinemet in y
#if AMREX_SPACEDIM == 3
                        int kk = iv[2] << 1; // refinemet in z
#endif
                        // iterate over all fine cells at the interface
                        // start with lower cells --> cover coarse neighbour
                        // cell
                        fill(map, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                        break;
                    }
                    case 1:
                    default:
                    {
                        // --> interface is on the upper face
                        int ii = covered[0] << 1; // refinemet in x
                        int jj = covered[1] << 1; // refinemet in y
#if AMREX_SPACEDIM == 3
                        int kk = covered[2] << 1; // refinemet in z
#endif
                        fill(map, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                        break;
                    }
                }
            }
        }
    }

    this->map2vector_m(map, indices, values);

    // iv: not covered coarse cell at crse-fine interface

    mglevel_m[level]->Bfine_p->insertGlobalValues(gidx,
                                                  indices.size(),
                                                  &values[0],
                                                  &indices[0]);
}


inline void AmrMultiGrid::buildDensityVector_m(const lo_t& level,
                                               const AmrField_t& rho)
{
    this->amrex2trilinos_m(level, 0, rho, mglevel_m[level]->rho_p);
}


inline void AmrMultiGrid::buildPotentialVector_m(const lo_t& level,
                                                 const AmrField_t& phi)
{
    this->amrex2trilinos_m(level, 0, phi, mglevel_m[level]->phi_p);
}


void AmrMultiGrid::buildGradientMatrix_m(const lo_t& level,
                                         const go_t& gidx,
                                         const AmrIntVect_t& iv,
                                         const basefab_t& mfab,
                                         const scalar_t* invdx)
{
    umap_t map;
    indices_t indices;
    coefficients_t values;

    auto check = [&](const AmrIntVect_t& iv,
                     const basefab_t& mfab,
                     int dir,
                     scalar_t shift)
    {
        switch ( mfab(iv) )
        {
            case AmrMultiGridLevel_t::Mask::INTERIOR:
                // interior cells
            case AmrMultiGridLevel_t::Mask::COVERED:
                // covered --> interior cell
                map[mglevel_m[level]->serialize(iv)] -= shift * 0.5 * invdx[dir];
                break;
            case AmrMultiGridLevel_t::Mask::BNDRY:
            {
                // interior boundary cells --> only level > 0
#ifndef NDEBUG
                if ( level == lbase_m )
                    throw OpalException("AmrMultiGrid::buildGradientMatrix_m()",
                                        "Error in mask for level "
                                        + std::to_string(level) + "!");
#endif

                scalar_t value = - shift * 0.5 * invdx[dir];

                // use 1st order Lagrange --> only cells of this level required
                AmrIntVect_t tmp = iv;
                // first fine cell on refined coarse cell (closer to interface)
                tmp[dir] -= shift;
                map[mglevel_m[level]->serialize(tmp)] += 2.0 * value;

                // second fine cell on refined coarse cell (further away from interface)
                tmp[dir] -= shift;
                map[mglevel_m[level]->serialize(tmp)] -= value;
                break;
            }
            case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
            {
                // physical boundary cells

                scalar_t value = - shift * 0.5 * invdx[dir];

                mglevel_m[level]->applyBoundary(iv, dir,
                                                map, value);
                break;
            }
            default:
                break;
        }
    };

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            AmrIntVect_t niv = iv;
            niv[d] += shift;
            check(niv, mfab, d, shift);
        }

        this->map2vector_m(map, indices, values);

        mglevel_m[level]->G_p[d]->insertGlobalValues(gidx,
                                                     indices.size(),
                                                     &values[0],
                                                     &indices[0]);
    }
}


void AmrMultiGrid::amrex2trilinos_m(const lo_t& level,
                                    const lo_t& comp,
                                    const AmrField_t& mf,
                                    Teuchos::RCP<vector_t>& mv)
{
    if ( mv.is_null() )
        mv = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );

    for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi) {
        const amrex::Box&          tbx  = mfi.tilebox();
        const amrex::FArrayBox&    fab = mf[mfi];

        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();

        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));

                    go_t gidx = mglevel_m[level]->serialize(iv);

                    mv->replaceGlobalValue(gidx, fab(iv, comp));
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
}


void AmrMultiGrid::trilinos2amrex_m(const lo_t& level,
                                    const lo_t& comp,
                                    AmrField_t& mf,
                                    const Teuchos::RCP<vector_t>& mv)
{
    Teuchos::ArrayRCP<const amr::scalar_t> data =  mv->get1dView();

    for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi) {
        const amrex::Box&          tbx  = mfi.tilebox();
        amrex::FArrayBox&          fab = mf[mfi];

        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();

        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));

                    go_t gidx = mglevel_m[level]->serialize(iv);
                    lo_t lidx = mglevel_m[level]->map_p->getLocalElement(gidx);

                    fab(iv, comp) = data[lidx];
                }
#if AMREX_SPACEDIM == 3
            }
#endif
        }
    }
}


inline
void AmrMultiGrid::map2vector_m(umap_t& map, indices_t& indices,
                                coefficients_t& values)
{
    indices.clear();
    values.clear();

    indices.reserve(map.size());
    values.reserve(map.size());

    std::for_each(map.begin(), map.end(),
                  [&](const std::pair<const go_t, scalar_t>& entry)
                  {
#ifndef NDEBUG
                      if ( entry.first < 0 ) {
                          throw OpalException("AmrMultiGrid::map2vector_m()",
                                              "Negative matrix index!");
                      }
#endif

                      indices.push_back(entry.first);
                      values.push_back(entry.second);
                  }
    );

    map.clear();
}


void AmrMultiGrid::smooth_m(const lo_t& level,
                            Teuchos::RCP<vector_t>& e,
                            Teuchos::RCP<vector_t>& r)
{
#if AMR_MG_TIMER
    IpplTimings::startTimer(smoothTimer_m);
#endif

    // base level has no smoother --> l - 1
    smoother_m[level-1]->smooth(e, mglevel_m[level]->Anf_p, r);

#if AMR_MG_TIMER
    IpplTimings::stopTimer(smoothTimer_m);
#endif
}


void AmrMultiGrid::restrict_m(const lo_t& level) {

#if AMR_MG_TIMER
        IpplTimings::startTimer(restrictTimer_m);
#endif

    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );

    this->residual_no_fine_m(level, tmp,
                             mglevel_m[level]->error_p,
                             mglevel_m[level-1]->error_p,
                             mglevel_m[level]->residual_p);

    mglevel_m[level-1]->residual_p->putScalar(0.0);

    // average down: residual^(l-1) = R^(l) * tmp
    mglevel_m[level-1]->R_p->apply(*tmp, *mglevel_m[level-1]->residual_p);


    // composite matrix, i.e. matrix without covered cells
    // r^(l-1) = rho^(l-1) - A * phi^(l-1)

    vector_t fine2crse(mglevel_m[level-1]->Awf_p->getDomainMap(), true);

    // get boundary for
    mglevel_m[level-1]->Bfine_p->apply(*mglevel_m[level]->phi_p, fine2crse);

    // operation: fine2coarse += A * phi
    mglevel_m[level-1]->Awf_p->apply(*mglevel_m[level-1]->phi_p,
                                     fine2crse, Teuchos::NO_TRANS,
                                     scalar_t(1.0), scalar_t(1.0));

    if ( mglevel_m[level-1]->Bcrse_p != Teuchos::null ) {
        //  operation: fine2coarse += B * phi
        mglevel_m[level-1]->Bcrse_p->apply(*mglevel_m[level-2]->phi_p,
                                           fine2crse, Teuchos::NO_TRANS,
                                           scalar_t(1.0), scalar_t(1.0));
    }

    Teuchos::RCP<vector_t> uncoveredRho = Teuchos::rcp( new vector_t(mglevel_m[level-1]->map_p, true) );

    mglevel_m[level-1]->UnCovered_p->apply(*mglevel_m[level-1]->rho_p, *uncoveredRho);


    //FIXME tmp2 not needed
    Teuchos::RCP<vector_t> tmp2 = Teuchos::rcp( new vector_t(mglevel_m[level-1]->map_p, true) );
    mglevel_m[level-1]->UnCovered_p->apply(fine2crse, *tmp2);

    // ONLY subtract coarse rho
    mglevel_m[level-1]->residual_p->update(1.0, *uncoveredRho, -1.0, *tmp2, 1.0);

#if AMR_MG_TIMER
        IpplTimings::stopTimer(restrictTimer_m);
#endif
}


void AmrMultiGrid::prolongate_m(const lo_t& level) {
#if AMR_MG_TIMER
        IpplTimings::startTimer(interpTimer_m);
#endif
        // interpolate error from l-1 to l
        // operation: e^(l) = 1.0 * e^(l) + 1.0 * I^(l) * e^(l-1)
        mglevel_m[level]->I_p->apply(*mglevel_m[level-1]->error_p,
                                     *mglevel_m[level]->error_p,
                                     Teuchos::NO_TRANS,
                                     scalar_t(1.0),
                                     scalar_t(1.0));
#if AMR_MG_TIMER
        IpplTimings::stopTimer(interpTimer_m);
#endif
}


void AmrMultiGrid::averageDown_m(const lo_t& level) {

    if (level == lfine_m )
        return;
#if AMR_MG_TIMER
        IpplTimings::startTimer(averageTimer_m);
#endif
    Teuchos::RCP<vector_t> phicrse = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );

    // operation: phicrse = 0.0 * phicrse + 1.0 * R^(l) * phi^(l+1)
    mglevel_m[level]->R_p->apply(*mglevel_m[level+1]->phi_p, *phicrse);

    Teuchos::RCP<vector_t> uncov_phi = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );

    mglevel_m[level]->UnCovered_p->apply(*mglevel_m[level]->phi_p, *uncov_phi);

    mglevel_m[level]->phi_p->update(1.0, *phicrse, 1.0, *uncov_phi, 0.0);
#if AMR_MG_TIMER
        IpplTimings::stopTimer(averageTimer_m);
#endif
}


void AmrMultiGrid::initInterpolater_m(const Interpolater& interp) {
    switch ( interp ) {
        case Interpolater::TRILINEAR:
            interp_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            throw OpalException("AmrMultiGrid::initInterpolater_m()",
                                "Not yet implemented.");
        case Interpolater::PIECEWISE_CONST:
            interp_mp.reset( new AmrPCInterpolater<AmrMultiGridLevel_t>() );
            break;
        default:
            throw OpalException("AmrMultiGrid::initInterpolater_m()",
                                "No such interpolater available.");
    }
}


void AmrMultiGrid::initCrseFineInterp_m(const Interpolater& interface) {
    switch ( interface ) {
        case Interpolater::TRILINEAR:
            interface_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            interface_mp.reset( new AmrLagrangeInterpolater<AmrMultiGridLevel_t>(
                AmrLagrangeInterpolater<AmrMultiGridLevel_t>::Order::QUADRATIC) );
            break;
        case Interpolater::PIECEWISE_CONST:
            interface_mp.reset( new AmrPCInterpolater<AmrMultiGridLevel_t>() );
            break;
        default:
            throw OpalException("AmrMultiGrid::initCrseFineInterp_m()",
                                "No such interpolater for the coarse-fine interface available.");
    }
}


void AmrMultiGrid::initBaseSolver_m(const BaseSolver& solver,
                                    const bool& rebalance,
                                    const std::string& reuse)
{
    switch ( solver ) {
        // Belos solvers
        case BaseSolver::BICGSTAB:
            solver_mp.reset( new BelosSolver_t("BICGSTAB", prec_mp) );
            break;
        case BaseSolver::MINRES:
            solver_mp.reset( new BelosSolver_t("MINRES", prec_mp) );
            break;
        case BaseSolver::PCPG:
            solver_mp.reset( new BelosSolver_t("PCPG", prec_mp) );
            break;
        case BaseSolver::CG:
            solver_mp.reset( new BelosSolver_t("Pseudoblock CG", prec_mp) );
            break;
        case BaseSolver::GMRES:
            solver_mp.reset( new BelosSolver_t("Pseudoblock GMRES", prec_mp) );
            break;
        case BaseSolver::STOCHASTIC_CG:
            solver_mp.reset( new BelosSolver_t("Stochastic CG", prec_mp) );
            break;
        case BaseSolver::RECYCLING_CG:
            solver_mp.reset( new BelosSolver_t("RCG", prec_mp) );
            break;
        case BaseSolver::RECYCLING_GMRES:
            solver_mp.reset( new BelosSolver_t("GCRODR", prec_mp) );
            break;
        // Amesos2 solvers
#ifdef HAVE_AMESOS2_KLU2
        case BaseSolver::KLU2:
            solver_mp.reset( new Amesos2Solver_t("klu2") );
            break;
#endif
#if HAVE_AMESOS2_SUPERLU
        case BaseSolver::SUPERLU:
            solver_mp.reset( new Amesos2Solver_t("superlu") );
            break;
#endif
#ifdef HAVE_AMESOS2_UMFPACK
        case BaseSolver::UMFPACK:
            solver_mp.reset( new Amesos2Solver_t("umfpack") );
            break;
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
        case BaseSolver::PARDISO_MKL:
            solver_mp.reset( new Amesos2Solver_t("pardiso_mkl") );
            break;
#endif
#ifdef HAVE_AMESOS2_MUMPS
        case BaseSolver::MUMPS:
            solver_mp.reset( new Amesos2Solver_t("mumps") );
            break;
#endif
#ifdef HAVE_AMESOS2_LAPACK
        case BaseSolver::LAPACK:
            solver_mp.reset( new Amesos2Solver_t("lapack") );
            break;
#endif
        case BaseSolver::SA:
        {
            std::string muelu = MueLuSolver_t::convertToMueLuReuseOption(reuse);
            solver_mp.reset( new MueLuSolver_t(rebalance, muelu) );
            break;
        }
        default:
            throw OpalException("AmrMultiGrid::initBaseSolver_m()",
                                "No such bottom solver available.");
    }
}


void AmrMultiGrid::initPrec_m(const Preconditioner& prec,
                              const bool& rebalance,
                              const std::string& reuse)
{
    switch ( prec ) {
        case Preconditioner::ILUT:
        case Preconditioner::CHEBYSHEV:
        case Preconditioner::RILUK:
        case Preconditioner::JACOBI:
        case Preconditioner::BLOCK_JACOBI:
        case Preconditioner::GS:
        case Preconditioner::BLOCK_GS:
            prec_mp.reset( new Ifpack2Preconditioner_t(prec) );
            break;
        case Preconditioner::SA:
        {
            std::string muelu = MueLuPreconditioner_t::convertToMueLuReuseOption(reuse);
            prec_mp.reset( new MueLuPreconditioner_t(rebalance, muelu) );
            break;
        }
        case Preconditioner::NONE:
            prec_mp.reset();
            break;
        default:
            throw OpalException("AmrMultiGrid::initPrec_m()",
                                "No such preconditioner available.");
    }
}


AmrMultiGrid::Boundary
AmrMultiGrid::convertToEnumBoundary_m(const std::string& bc) {
    std::map<std::string, Boundary> map;

    map["DIRICHLET"]    = Boundary::DIRICHLET;
    map["OPEN"]         = Boundary::OPEN;
    map["PERIODIC"]     = Boundary::PERIODIC;

    auto boundary = map.find(Util::toUpper(bc));

    if ( boundary == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumBoundary_m()",
                            "No boundary type '" + bc + "'.");
    return boundary->second;
}

AmrMultiGrid::Interpolater
AmrMultiGrid::convertToEnumInterpolater_m(const std::string& interp) {
    std::map<std::string, Interpolater> map;

    map["TRILINEAR"]    = Interpolater::TRILINEAR;
    map["LAGRANGE"]     = Interpolater::LAGRANGE;
    map["PC"]           = Interpolater::PIECEWISE_CONST;

    auto interpolater = map.find(Util::toUpper(interp));

    if ( interpolater == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumInterpolater_m()",
                            "No interpolater '" + interp + "'.");
    return interpolater->second;
}


AmrMultiGrid::BaseSolver
AmrMultiGrid::convertToEnumBaseSolver_m(const std::string& bsolver) {
    std::map<std::string, BaseSolver> map;

    map["BICGSTAB"]         = BaseSolver::BICGSTAB;
    map["MINRES"]           = BaseSolver::MINRES;
    map["PCPG"]             = BaseSolver::PCPG;
    map["CG"]               = BaseSolver::CG;
    map["GMRES"]            = BaseSolver::GMRES;
    map["STOCHASTIC_CG"]    = BaseSolver::STOCHASTIC_CG;
    map["RECYCLING_CG"]     = BaseSolver::RECYCLING_GMRES;
    map["RECYCLING_GMRES"]  = BaseSolver::RECYCLING_GMRES;
#ifdef HAVE_AMESOS2_KLU2
    map["KLU2"]             = BaseSolver::KLU2;
#endif
#ifdef HAVE_AMESOS2_SUPERLU
    map["SUPERLU"]          = BaseSolver::SUPERLU;
#endif
#ifdef HAVE_AMESOS2_UMFPACK
    map["UMFPACK"]          = BaseSolver::UMFPACK;
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
    map["PARDISO_MKL"]      = BaseSolver::PARDISO_MKL;
#endif
#ifdef HAVE_AMESOS2_MUMPS
    map["MUMPS"]            = BaseSolver::MUMPS;
#endif
#ifdef HAVE_AMESOS2_LAPACK
    map["LAPACK"]           = BaseSolver::LAPACK;
#endif
    map["SA"]               = BaseSolver::SA;

    auto solver = map.find(Util::toUpper(bsolver));

    if ( solver == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumBaseSolver_m()",
                            "No bottom solver '" + bsolver + "'.");
    return solver->second;
}


AmrMultiGrid::Preconditioner
AmrMultiGrid::convertToEnumPreconditioner_m(const std::string& prec) {
    std::map<std::string, Preconditioner> map;

    map["NONE"]         = Preconditioner::NONE;

    Ifpack2Preconditioner_t::fillMap(map);

    MueLuPreconditioner_t::fillMap(map);

    auto precond = map.find(Util::toUpper(prec));

    if ( precond == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumPreconditioner_m()",
                            "No preconditioner '" + prec + "'.");
    return precond->second;
}


AmrMultiGrid::Smoother
AmrMultiGrid::convertToEnumSmoother_m(const std::string& smoother) {
    return AmrSmoother::convertToEnumSmoother(smoother);
}


AmrMultiGrid::Norm
AmrMultiGrid::convertToEnumNorm_m(const std::string& norm) {
    std::map<std::string, Norm> map;

    map["L1"]   = Norm::L1;
    map["L2"]   = Norm::L2;
    map["LINF"] = Norm::LINF;

    snorm_m = Util::toUpper(norm);

    auto n = map.find(snorm_m);

    if ( n == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumNorm_m()",
                            "No norm '" + norm + "'.");
    return n->second;
}


void AmrMultiGrid::writeSDDSHeader_m(std::ofstream& outfile) {
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    outfile << "SDDS1" << std::endl;
    outfile << "&description\n"
            << indent << "text=\"Solver statistics '" << OpalData::getInstance()->getInputFn()
            << "' " << dateStr << "" << timeStr << "\",\n"
            << indent << "contents=\"solver info\"\n"
            << "&end\n";
    outfile << "&parameter\n"
            << indent << "name=processors,\n"
            << indent << "type=long,\n"
            << indent << "description=\"Number of Cores used\"\n"
            << "&end\n";
    outfile << "&parameter\n"
            << indent << "name=revision,\n"
            << indent << "type=string,\n"
            << indent << "description=\"git revision of opal\"\n"
            << "&end\n";
    outfile << "&parameter\n"
            << indent << "name=flavor,\n"
            << indent << "type=string,\n"
            << indent << "description=\"OPAL flavor that wrote file\"\n"
            << "&end\n";
    outfile << "&column\n"
            << indent << "name=t,\n"
            << indent << "type=double,\n"
            << indent << "units=ns,\n"
            << indent << "description=\"1 Time\"\n"
            << "&end\n";
    outfile << "&column\n"
            << indent << "name=mg_iter,\n"
            << indent << "type=long,\n"
            << indent << "units=1,\n"
            << indent << "description=\"2 Number of Multigrid Iterations\"\n"
            << "&end\n";
    outfile << "&column\n"
            << indent << "name=bottom_iter,\n"
            << indent << "type=long,\n"
            << indent << "units=1,\n"
            << indent << "description=\"3 Total Number of Bottom Solver Iterations\"\n"
            << "&end\n";
    outfile << "&column\n"
            << indent << "name=regrid,\n"
            << indent << "type=bool,\n"
            << indent << "units=1,\n"
            << indent << "description=\"4 Regrid Step\"\n"
            << "&end\n";
    outfile << "&column\n"
            << indent << "name=" + snorm_m + ",\n"
            << indent << "type=double,\n"
            << indent << "units=1,\n"
            << indent << "description=\"5 Error\"\n"
            << "&end\n"
            << "&data\n"
            << indent << "mode=ascii,\n"
            << indent << "no_row_counts=1\n"
            << "&end\n"
            << comm_mp->getSize() << '\n'
            << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << '\n'
            << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void AmrMultiGrid::writeSDDSData_m(const scalar_t& error) {
#if AMR_MG_TIMER
    IpplTimings::startTimer(dumpTimer_m);
#endif
    unsigned int pwi = 10;

    std::ofstream outfile;

    if ( comm_mp->getRank() == 0 ) {
        outfile.open(fname_m.c_str(), flag_m);
        outfile.precision(15);
        outfile.setf(std::ios::scientific, std::ios::floatfield);

        if ( flag_m == std::ios::out ) {
            flag_m = std::ios::app;
            writeSDDSHeader_m(outfile);
        }

        outfile << itsAmrObject_mp->getT() * 1e9 << std::setw(pwi) << '\t'  // 1
                << this->nIter_m << std::setw(pwi) << '\t'                  // 2
                << this->bIter_m << std::setw(pwi) << '\t'                  // 3
                << this->regrid_m << std::setw(pwi) <<  '\t'                // 4
                << error << '\n';                                           // 5
    }
#if AMR_MG_TIMER
    IpplTimings::stopTimer(dumpTimer_m);
#endif
}


#if AMR_MG_TIMER
void AmrMultiGrid::initTimer_m() {
    buildTimer_m        = IpplTimings::getTimer("AMR MG matrix setup");
    restrictTimer_m     = IpplTimings::getTimer("AMR MG restrict");
    smoothTimer_m       = IpplTimings::getTimer("AMR MG smooth");
    interpTimer_m       = IpplTimings::getTimer("AMR MG prolongate");
    efieldTimer_m       = IpplTimings::getTimer("AMR MG e-field");
    averageTimer_m      = IpplTimings::getTimer("AMR MG average down");
    bottomTimer_m       = IpplTimings::getTimer("AMR MG bottom-solver");
    dumpTimer_m         = IpplTimings::getTimer("AMR MG dump");
}
#endif


double AmrMultiGrid::getXRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(0);
}


double AmrMultiGrid::getXRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(0);
}


double AmrMultiGrid::getYRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(1);
}


double AmrMultiGrid::getYRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(1);
}


double AmrMultiGrid::getZRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(2);
}


double AmrMultiGrid::getZRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(2);
}


Inform &AmrMultiGrid::print(Inform &os) const {
    os << "* ********************* A M R M u l t i G r i d ********************** " << endl
    //FIXME
       << "* ******************************************************************** " << endl;
    return os;
}