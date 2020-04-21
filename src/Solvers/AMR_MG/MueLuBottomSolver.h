//
// Class MueLuBottomSolver
//   Interface to the SAAMG solver of MueLu.
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
#ifndef MUELU_BOTTOM_SOLVER_H
#define MUELU_BOTTOM_SOLVER_H

#include "BottomSolver.h"

#include "Amr/AmrDefs.h"

#include "Ippl.h"

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Utilities.hpp>

#include <MueLu_ParameterListInterpreter.hpp>

template <class Level>
class MueLuBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                              Teuchos::RCP<amr::multivector_t>,
                                              Level >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::node_t node_t;

    typedef amr::AmrBox_t AmrBox_t;
    typedef amr::AmrIntVect_t AmrIntVect_t;

//    typedef amr::AmrGeometry_t AmrGeometry_t;
    
    typedef MueLu::Hierarchy<scalar_t, lo_t, go_t, node_t> hierarchy_t;
    typedef MueLu::Level level_t;
    typedef Xpetra::Matrix<scalar_t, lo_t, go_t, node_t> xmatrix_t;
    typedef Xpetra::MultiVector<scalar_t, lo_t, go_t, node_t> xmv_t;
    typedef MueLu::Utilities<scalar_t, lo_t, go_t, node_t> util_t;
    
    typedef MueLu::ParameterListInterpreter<scalar_t, lo_t, go_t, node_t> pListInterpreter_t;
    typedef MueLu::HierarchyManager<scalar_t, lo_t, go_t, node_t> manager_t;
        
public:
    
    MueLuBottomSolver(const bool& rebalance,
                      const std::string& reuse);
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A,
                     Level* level_p = nullptr);
    
    std::size_t getNumIters();
    
    /*
     * MueLu reuse option.
     * Either none, RP, RAP or full
     */
    static std::string convertToMueLuReuseOption(const std::string& reuse);

private:
    void initMueLuList_m(const std::string& reuse);
    
private:
    Teuchos::RCP<hierarchy_t> hierarchy_mp;     ///< manages the multigrid hierarchy
    
    Teuchos::RCP<manager_t> factory_mp;         ///< sets up hierarchy

    Teuchos::RCP<xmatrix_t> A_mp;               ///< MueLu requires Xpetra

    lo_t nSweeps_m;                             ///< the number of multigrid iterations
    
    Teuchos::ParameterList mueluList_m;
    
    bool rebalance_m;                           ///< use subcommunicators (less communication)
    
    IpplTimings::TimerRef setupTimer_m;
};

#include "MueLuBottomSolver.hpp"

#endif
