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
#ifndef MUELU_PRECONDITIONER_H
#define MUELU_PRECONDITIONER_H

#include "AmrPreconditioner.h"
#include "Amr/AmrDefs.h"

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>

template <class Level>
class MueLuPreconditioner : public AmrPreconditioner<amr::matrix_t, Level>
{
public:
    typedef amr::Preconditioner Preconditioner;
    
    typedef amr::scalar_t scalar_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::AmrBox_t AmrBox_t;
    
    typedef MueLu::TpetraOperator<
        scalar_t,
        lo_t,
        go_t,
        amr::node_t
    > precond_t;
    
    typedef amr::AmrIntVect_t AmrIntVect_t;
    
    typedef std::map<std::string, Preconditioner> map_t;
    
public:
    
    MueLuPreconditioner(const bool& rebalance,
                        const std::string& reuse);
    
    void create(const Teuchos::RCP<amr::matrix_t>& A, Level* level_p =  nullptr);
    
    Teuchos::RCP<amr::operator_t> get();
    
    static void fillMap(map_t& map);
    
    static std::string convertToMueLuReuseOption(const std::string& reuse);
    
private:
    void init_m(const std::string& reuse);

private:
    Teuchos::ParameterList params_m;
    
    Teuchos::RCP<precond_t> prec_mp;

    Teuchos::RCP<amr::multivector_t> coords_mp;
    
    const bool rebalance_m;
};

#include "MueLuPreconditioner.hpp"

#endif
