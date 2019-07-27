/* https://trilinos.org/wordpress/wp-content/uploads/2015/03/MueLu_tutorial.pdf
 * http://prod.sandia.gov/techlib/access-control.cgi/2014/1418624r.pdf
 */

#include <AMReX.H>
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

template <class Level>
MueLuBottomSolver<Level>::MueLuBottomSolver(const bool& rebalance,
                                            const std::string& reuse)
    : BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                Teuchos::RCP<amr::multivector_t>,
                                Level>()
    , A_mp(Teuchos::null)
    , nSweeps_m(4)
    , rebalance_m(rebalance)
    , setupTimer_m(IpplTimings::getTimer("AMR MG bsolver setup"))
{
    this->initMueLuList_m(reuse);
    
    factory_mp = Teuchos::rcp( new pListInterpreter_t(mueluList_m) );
    
    // empty multigrid hierarchy with a finest level only
    hierarchy_mp = factory_mp->CreateHierarchy();
}


template <class Level>
void MueLuBottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                                     const Teuchos::RCP<mv_t>& b)
{
    // MueLu requires Xpetra multivectors (wrap them)
    Teuchos::RCP<xmv_t> xx = MueLu::TpetraMultiVector_To_XpetraMultiVector(x);
    Teuchos::RCP<xmv_t> xb = MueLu::TpetraMultiVector_To_XpetraMultiVector(b);

    // InitialGuessIsZero = true
    // startLevel = 0
    hierarchy_mp->Iterate(*xb, *xx, nSweeps_m, true, 0);
    
    // put multivector back
    x->assign(*util_t::MV2NonConstTpetraMV2(*xx));
}


template <class Level>
void MueLuBottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A,
                                           Level* level_p)
{
    IpplTimings::startTimer(setupTimer_m);

    A_mp = MueLu::TpetraCrs_To_XpetraMatrix<scalar_t, lo_t, go_t, node_t>(A);
    A_mp->SetFixedBlockSize(1); // only 1 DOF per node (pure Laplace problem)

    Teuchos::RCP<mv_t> coords_p = Teuchos::rcp(
        new amr::multivector_t(A->getDomainMap(), AMREX_SPACEDIM, false)
    );

    const scalar_t* domain = level_p->geom.ProbLo();
    const scalar_t* dx = level_p->cellSize();
    for (amrex::MFIter mfi(level_p->grids, level_p->dmap, true);
        mfi.isValid(); ++mfi)
    {
        const AmrBox_t&       tbx = mfi.tilebox();
        const lo_t* lo = tbx.loVect();
        const lo_t* hi = tbx.hiVect();

        for (lo_t i = lo[0]; i <= hi[0]; ++i) {
            for (lo_t j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (lo_t k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    go_t gidx = level_p->serialize(iv);

                    coords_p->replaceGlobalValue(gidx, 0, domain[0] + (0.5 + i) * dx[0]);
                    coords_p->replaceGlobalValue(gidx, 1, domain[1] + (0.5 + j) * dx[1]);
#if AMREX_SPACEDIM == 3
                    coords_p->replaceGlobalValue(gidx, 2, domain[2] + (0.5 + k) * dx[2]);
                }
#endif
            }
        }
    }

    Teuchos::RCP<xmv_t> coordinates = MueLu::TpetraMultiVector_To_XpetraMultiVector(coords_p);

    Teuchos::RCP<mv_t> nullspace = Teuchos::rcp(new mv_t(A->getRowMap(), 1));
    Teuchos::RCP<xmv_t> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector(nullspace);
    xnullspace->putScalar(1.0);
    hierarchy_mp->GetLevel(0)->Set("Nullspace", xnullspace);
    hierarchy_mp->GetLevel(0)->Set("Coordinates", coordinates);
    hierarchy_mp->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    hierarchy_mp->IsPreconditioner(false);
    hierarchy_mp->GetLevel(0)->Set("A", A_mp);


    Teuchos::RCP<level_t> finest_p = hierarchy_mp->GetLevel(0);
    finest_p->Set("A", A_mp);

    factory_mp->SetupHierarchy(*hierarchy_mp);

    this->isInitialized_m = true;

    IpplTimings::stopTimer(setupTimer_m);
}


template <class Level>
std::size_t MueLuBottomSolver<Level>::getNumIters() {
    return nSweeps_m;
}


template <class Level>
std::string
MueLuBottomSolver<Level>::convertToMueLuReuseOption(const std::string& reuse) {
    
    std::map<std::string, std::string> map;
    map["NONE"] = "none";
    map["RP"]   = "RP";
    map["RAP"]  = "RAP";
    map["S"]    = "S";
    map["FULL"] = "full";
    
    auto muelu =  map.find(Util::toUpper(reuse));
    
    if ( muelu == map.end() )
        throw OpalException("MueLuBottomSolver::convertToMueLuReuseOption()",
                            "No MueLu reuse option '" + reuse + "'.");
    
    return muelu->second;
}


template <class Level>
void MueLuBottomSolver<Level>::initMueLuList_m(const std::string& reuse) {
    mueluList_m.set("problem: type", "Poisson-3D");
    mueluList_m.set("verbosity", "low");
    mueluList_m.set("number of equations", 1);
    mueluList_m.set("max levels", 8);
    mueluList_m.set("cycle type", "V");

    mueluList_m.set("coarse: max size", 200);
    mueluList_m.set("multigrid algorithm", "sa");
    mueluList_m.set("sa: damping factor", 1.33); // default: 1.33
    mueluList_m.set("sa: use filtered matrix", true);
    mueluList_m.set("filtered matrix: reuse eigenvalue", false); // false: more expensive
    
    mueluList_m.set("repartition: enable", rebalance_m);
    mueluList_m.set("repartition: rebalance P and R", rebalance_m);
    mueluList_m.set("repartition: partitioner", "zoltan2");
    mueluList_m.set("repartition: min rows per proc", 800);
    mueluList_m.set("repartition: start level", 2);

    Teuchos::ParameterList reparms;
    reparms.set("algorithm", "multijagged"); // rcb
    //    reparms.set("partitioning_approach", "partition");

    mueluList_m.set("repartition: params", reparms);
    
    mueluList_m.set("smoother: type", "CHEBYSHEV");
    mueluList_m.set("smoother: pre or post", "both");
    Teuchos::ParameterList smparms;
    smparms.set("chebyshev: degree", 3);
    smparms.set("chebyshev: assume matrix does not change", false);
    smparms.set("chebyshev: zero starting solution", true);
    mueluList_m.set("smoother: params", smparms);


    mueluList_m.set("coarse: type", "RELAXATION");
    Teuchos::ParameterList cparms;
    cparms.set("relaxation: type", "Gauss-Seidel");
    cparms.set("relaxation: sweeps", nSweeps_m);
    cparms.set("relaxation: zero starting solution", true);
    cparms.set("relaxation: use l1", true);
    cparms.set("relaxation: l1 eta", 1.5);
    mueluList_m.set("coarse: params", cparms);

//    mueluList_m.set("coarse: type", "KLU2");

    mueluList_m.set("aggregation: type", "uncoupled");
    mueluList_m.set("aggregation: min agg size", 3);
    mueluList_m.set("aggregation: max agg size", 27);

    mueluList_m.set("transpose: use implicit", false);

    mueluList_m.set("reuse: type", reuse);
}

