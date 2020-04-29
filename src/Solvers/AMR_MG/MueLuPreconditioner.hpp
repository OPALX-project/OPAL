//
// Class MueLuPreconditioner
//   Interface to the SAAMG solver of MueLu. Here it is used as preconditioner for Belos
//   iterative solvers.
//
// Copyright (c) 2017 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include <AMReX.H>

#include <MueLu_CreateTpetraPreconditioner.hpp>

template <class Level>
MueLuPreconditioner<Level>::MueLuPreconditioner(const bool& rebalance,
                                                const std::string& reuse)
    : prec_mp(Teuchos::null),
      coords_mp(Teuchos::null),
      rebalance_m(rebalance)
{
    this->init_m(reuse);
}


template <class Level>
void MueLuPreconditioner<Level>::create(const Teuchos::RCP<amr::matrix_t>& A,
                                        Level* level_p)
{

    coords_mp = Teuchos::null;

    if ( rebalance_m ) {
        
        const scalar_t* domain = level_p->geom.ProbLo();
        const scalar_t* dx = level_p->cellSize();
        
        coords_mp = Teuchos::rcp(
            new amr::multivector_t(A->getDomainMap(), AMREX_SPACEDIM, false)
        );
    
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
                    
                        coords_mp->replaceGlobalValue(gidx, 0, domain[0] + (0.5 + i) * dx[0]);
                        coords_mp->replaceGlobalValue(gidx, 1, domain[1] + (0.5 + j) * dx[1]);
#if AMREX_SPACEDIM == 3
                        coords_mp->replaceGlobalValue(gidx, 2, domain[2] + (0.5 + k) * dx[2]);
                    }
#endif
                }
            }
        }
    }

    prec_mp = MueLu::CreateTpetraPreconditioner(A, params_m, coords_mp);

    coords_mp = Teuchos::null;
}


template <class Level>
Teuchos::RCP<amr::operator_t> MueLuPreconditioner<Level>::get() {
    return prec_mp;
}


template <class Level>
void MueLuPreconditioner<Level>::fillMap(map_t& map) {
    map["SA"] = Preconditioner::SA;
}


template <class Level>
std::string
MueLuPreconditioner<Level>::convertToMueLuReuseOption(const std::string& reuse) {
    
    std::map<std::string, std::string> map;
    map["NONE"] = "none";
    map["RP"]   = "RP";
    map["RAP"]  = "RAP";
    map["S"]    = "S";
    map["FULL"] = "full";
    
    auto muelu =  map.find(reuse);
    
    if ( muelu == map.end() )
        throw OpalException("MueLuPreconditioner::convertToMueLuReuseOption()",
                            "No MueLu reuse option '" + reuse + "'.");
    
    return muelu->second;
}


template <class Level>
void MueLuPreconditioner<Level>::init_m(const std::string& reuse) {
    params_m.set("problem: type", "Poisson-3D");
    params_m.set("verbosity", "none");
    params_m.set("number of equations", 1);
    params_m.set("max levels", 8);
    params_m.set("cycle type", "V");

    params_m.set("coarse: max size", 200);
    params_m.set("multigrid algorithm", "sa");
    params_m.set("sa: damping factor", 1.33); // default: 1.33
    params_m.set("sa: use filtered matrix", true);
    params_m.set("filtered matrix: reuse eigenvalue", false); // false: more expensive
    
    params_m.set("repartition: enable", rebalance_m);
    params_m.set("repartition: rebalance P and R", rebalance_m);
    params_m.set("repartition: partitioner", "zoltan2");
    params_m.set("repartition: min rows per proc", 800);
    params_m.set("repartition: start level", 2);

    Teuchos::ParameterList reparms;
    reparms.set("algorithm", "multijagged"); // rcb
    //    reparms.set("partitioning_approach", "partition");

    params_m.set("repartition: params", reparms);
    
    
    params_m.set("smoother: type", "CHEBYSHEV");
    params_m.set("smoother: pre or post", "both");
    Teuchos::ParameterList smparms;
    smparms.set("chebyshev: degree", 3);
    smparms.set("chebyshev: assume matrix does not change", false);
    smparms.set("chebyshev: zero starting solution", true);
    params_m.set("smoother: params", smparms);
    
    params_m.set("smoother: type", "CHEBYSHEV");
    params_m.set("smoother: pre or post", "both");

    params_m.set("coarse: type", "KLU2");

    params_m.set("aggregation: type", "uncoupled");
    params_m.set("aggregation: min agg size", 3);
    params_m.set("aggregation: max agg size", 27);

    params_m.set("transpose: use implicit", false);

    params_m.set("reuse: type", reuse); // none
}
