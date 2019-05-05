#include "FMGPoissonSolver.h"

#include "Utilities/OpalException.h"

#include <AMReX_ParmParse.H>
#include <AMReX_Interpolater.H>

FMGPoissonSolver::FMGPoissonSolver(AmrBoxLib* itsAmrObject_p)
    : AmrPoissonSolver<AmrBoxLib>(itsAmrObject_p),
      reltol_m(1.0e-15),
      abstol_m(1.0e-9)
{
    // Dirichlet boundary conditions are default
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc_m[2 * d]     = MGT_BC_DIR;
        bc_m[2 * d + 1] = MGT_BC_DIR;
    }

    this->initParameters_m();
}

void FMGPoissonSolver::solve(AmrScalarFieldContainer_t& rho,
                             AmrScalarFieldContainer_t& phi,
                             AmrVectorFieldContainer_t& efield,
                             unsigned short baseLevel,
                             unsigned short finestLevel,
                             bool prevAsGuess)
{
    const GeomContainer_t& geom = itsAmrObject_mp->Geom();

    if (AmrGeometry_t::isAllPeriodic()) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bc_m[2 * d]     = MGT_BC_PER;
            bc_m[2 * d + 1] = MGT_BC_PER;
        }
    } else if ( AmrGeometry_t::isAnyPeriodic() ) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if ( AmrGeometry_t::isPeriodic(d) ) {
                bc_m[2 * d]     = MGT_BC_PER;
                bc_m[2 * d + 1] = MGT_BC_PER;
            }
        }
    }

    amrex::Vector< AmrScalarFieldContainer_t > grad_phi_edge(rho.size());

    for (int lev = baseLevel; lev <= finestLevel ; ++lev) {
        const AmrProcMap_t& dmap = rho[lev]->DistributionMap();
        grad_phi_edge[lev].resize(AMREX_SPACEDIM);

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            AmrGrid_t ba = rho[lev]->boxArray();
            grad_phi_edge[lev][n].reset(new AmrField_t(ba.surroundingNodes(n), dmap, 1, 1));
            grad_phi_edge[lev][n]->setVal(0.0, 1);
        }
    }

    for (int i = 0; i <= finestLevel; ++i) {
        phi[i]->setVal(0.0, 1);
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield[i][j]->setVal(0.0, 1);
        }
    }

    double residNorm = this->solveWithF90_m(amrex::GetVecOfPtrs(rho),
                                            amrex::GetVecOfPtrs(phi),
                                            amrex::GetVecOfVecOfPtrs(grad_phi_edge),
                                            geom,
                                            baseLevel,
                                            finestLevel);

    if ( residNorm > abstol_m ) {
        std::stringstream ss;
        ss << "Residual norm: " << std::setprecision(16) << residNorm
           << " > " << abstol_m << " (absolute tolerance)";
        throw OpalException("FMGPoissonSolver::solve()",
                            "Multigrid solver did not converge. " + ss.str());
    }

    for (int lev = baseLevel; lev <= finestLevel; ++lev) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            amrex::average_face_to_cellcenter(*(efield[lev][j].get()),
                                              amrex::GetVecOfConstPtrs(grad_phi_edge[lev]),
                                              geom[lev]);

            efield[lev][j]->FillBoundary(0, 1, geom[lev].periodicity());
            // we need also minus sign due to \vec{E} = - \nabla\phi
            efield[lev][j]->mult(-1.0, 0, 1);
        }
    }
}


double FMGPoissonSolver::getXRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(0);
}


double FMGPoissonSolver::getXRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(0);
}


double FMGPoissonSolver::getYRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(1);
}


double FMGPoissonSolver::getYRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(1);
}


double FMGPoissonSolver::getZRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(2);
}


double FMGPoissonSolver::getZRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(2);
}


Inform &FMGPoissonSolver::print(Inform &os) const {
    os << "* ************* F M G P o i s s o n S o l v e r ************************************ " << endl
       << "* relative tolerance " << reltol_m << '\n'
       << "* absolute tolerance " << abstol_m << '\n' << endl
       << "* ******************************************************************** " << endl;
    return os;
}


void FMGPoissonSolver::initParameters_m() {
    /* The information about the parameters is copied from the
     * BoxLib documentation chapter 5 on linear solvers and from
     * the code itself.
     *
     * Some paramters that have to be given to the solver using
     * AMReX_ParmParse.
     * "doc" tells how the variable is called in
     *  - amrex/Src/LinearSolvers/F_MG/cc_mg_cpp.f90
     *  - amrex/Src/LinearSolvers/F_MG/mg_tower.f90 (conains defaults)
     *  - amrex/Src/LinearSolvers/C_to_F_MG/AMReX_MGT_Solver.cpp
     *
     */

    amrex::ParmParse pp_mg("mg");

    // maximum number of multigrid cycles (doc: max_iter)
    pp_mg.add("maxiter", 200);

    // maximum number of iterations for the bottom solver (doc: bottom_max_iter)
    pp_mg.add("maxiter_b", 200);

    // number of smoothings at each level on the way DOWN the V-cycle (doc: nu1)
    pp_mg.add("nu_1", 2);

    // number of smoothings at each level on the way UP the V-cycle (doc: nu2)
    pp_mg.add("nu_2", 2);

    // number of smoothing before and after the bottom solver (doc: nub)
    pp_mg.add("nu_b", 0);

    // number of smoothing ... ?
    pp_mg.add("nu_f", 8);

    // verbosity of the multigrid solver. Higher numbers give more verbosity (doc: verbose)
    pp_mg.add("v"   , 0);


    // see amrex/Src/LinearSolvers/C_to_F_MG/AMReX_MGT_Solver.cpp
    pp_mg.add("usecg", 1);

    // bottom_solver_eps, see amrex/Src/LinearSolvers/C_to_F_MG/AMReX_MGT_Solver.cpp
    pp_mg.add("rtol_b", 0.0001);

    // see amrex/Src/LinearSolvers/C_to_F_MG/AMReX_MGT_Solver.cpp
    pp_mg.add("cg_solver", 1);

    // maximum number of levels (max_nlevel)
    pp_mg.add("numLevelsMAX", 1024);

    /* MG_SMOOTHER_GS_RB  = 1   (red-black Gauss-Seidel)
     * MG_SMOOTHER_JACOBI = 2   (Jacobi)
     * MG_SMOOTHER_MINION_CROSS = 5
     * MG_SMOOTHER_MINION_FULL = 6
     * MG_SMOOTHER_EFF_RB = 7
     */
    pp_mg.add("smoother", 1);

    /* The type of multigrid to do
     *
     * MG_FCycle = 1    (F, full multigrid)
     * MG_WCycle = 2    (W-cycles)
     * MG_VCycle = 3    (V-cycles)
     * MG_FVCycle = 4   (F + V cycles)
     */
    pp_mg.add("cycle_type", 1);

    /* the type of bottom solver to use
     *
     * BiCG:    biconjugate gradient stabilized (bottom_solver = 1)
     * CG:      conjugate gradient method (bottom_solver = 2)
     * CABiCG:  (bottom_solver = 3)
     * special: bottom_solver = 4
     *
     * The special bottom solver extends the range of the multigrid coarsening
     * by aggregating coarse grids on the original mesh together and further coarsening.
     *
     * if use_cg == 1 && cg_solver == 0 --> CG
     * if use_cg == 1 && cg_solver == 1 --> BiCG
     * if use_cg == 1 && cg_solver == 2 --> CABiCG
     * if use_cg == 0                   --> CABiCG
     */
    pp_mg.add("bottom_solver", 1);

    // verbosity of the bottom solver. Higher numbers give more verbosity (doc: cg_verbose)
    amrex::ParmParse pp_cg("cg");
    pp_cg.add("v", 0);


    /* Additional parameters that can't be set by ParmParse:
     *
     *      - max_bottom_nlevel = 3 (additional coarsening if you use bottom_solver == 4)
     *      - min_width = 2 (minimum size of grid at coarsest multigrid level)
     *      - max_L0_growth = -1
     */
    amrex::MGT_Solver::def_min_width = 2;
    amrex::MGT_Solver::def_max_L0_growth = -1.0;
}


double FMGPoissonSolver::solveWithF90_m(const AmrFieldContainer_pt& rho,
                                        const AmrFieldContainer_pt& phi,
                                        const amrex::Vector< AmrFieldContainer_pt > & grad_phi_edge,
                                        const GeomContainer_t& geom,
                                        int baseLevel,
                                        int finestLevel)
{
    int nlevs = finestLevel - baseLevel + 1;

    GeomContainer_t geom_p(nlevs);
    AmrFieldContainer_pt rho_p(nlevs);
    AmrFieldContainer_pt phi_p(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p[ilev] = geom[ilev + baseLevel];
        rho_p[ilev]  = rho[ilev + baseLevel];
        phi_p[ilev]  = phi[ilev + baseLevel];
    }

    //FIXME Refinement ratio
    amrex::IntVect crse_ratio = (baseLevel == 0) ?
        amrex::IntVect::TheZeroVector() : itsAmrObject_mp->refRatio(0);

    amrex::FMultiGrid fmg(geom_p, baseLevel, crse_ratio);

    if (baseLevel == 0)
        fmg.set_bc(bc_m, *phi_p[baseLevel]);
    else
        fmg.set_bc(bc_m, *phi_p[baseLevel-1], *phi_p[baseLevel]);

    /* (alpha * a - beta * (del dot b grad)) phi = rho
     * (b = 1)
     *
     * The function call set_const_gravity_coeffs() sets alpha = 0.0
     * and beta = -1 (in MGT_Solver::set_const_gravity_coeffs)
     *
     * --> (del dot grad) phi = rho
     */
    fmg.set_const_gravity_coeffs();

    // order of approximation at Dirichlet boundaries
    fmg.set_maxorder(3);

    int always_use_bnorm = 0;
    int need_grad_phi = 1;

    double residNorm = fmg.solve(phi_p, rho_p, reltol_m, abstol_m, always_use_bnorm, need_grad_phi);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + baseLevel;
        fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
    }

    return residNorm;
}


/*
void FMGPoissonSolver::interpolate_m(AmrScalarFieldContainer_t& phi,
                                     const GeomContainer_t& geom,
                                     double l0norm,
                                     int finestLevel)
{
    amrex::PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Vector<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;

    std::unique_ptr<AmrField_t> tmp;

    for (int lev = 1; lev <= finestLevel; ++lev) {
        const AmrGrid_t& ba = phi[lev]->boxArray();
        const AmrProcMap_t& dm = phi[lev]->DistributionMap();
        tmp.reset(new AmrField_t(ba, dm, 1, 0));
        tmp->setVal(0.0);

        amrex::InterpFromCoarseLevel(*tmp, 0.0, *phi[lev-1],
                                     0, 0, 1, geom[lev-1], geom[lev],
                                     cphysbc, fphysbc,
                                     itsAmrObject_mp->refRatio(lev-1),
                                     &mapper, bcs);
        phi[lev]->plus(*tmp, 0, 1, 0);
        phi[lev-1]->mult(l0norm, 0, 1);
    }

    phi[finestLevel]->mult(l0norm, 0, 1);
}
*/