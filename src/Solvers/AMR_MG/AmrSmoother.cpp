//
// Class AmrSmoother
//   Interface to Ifpack2 smoothers of the Trilinos package.
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
#include "AmrSmoother.h"

#include <map>
#include <utility>

#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

AmrSmoother::AmrSmoother(const Teuchos::RCP<const matrix_t>& A,
                         const Smoother& smoother,
                         lo_t nSweeps)
{
    const std::string type = "RELAXATION";
    
    Ifpack2::Factory factory;
    prec_mp = factory.create(type, A);
    
    params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    
    this->initParameter_m(smoother, nSweeps);
    
    
    prec_mp->setParameters(*params_mp);
    prec_mp->initialize();
    prec_mp->compute();
}


AmrSmoother::~AmrSmoother() {
    prec_mp = Teuchos::null;
    params_mp = Teuchos::null;
}


void AmrSmoother::smooth(const Teuchos::RCP<vector_t>& x,
                         const Teuchos::RCP<matrix_t>& A,
                         const Teuchos::RCP<vector_t>& b)
{
    prec_mp->apply(*b, *x, Teuchos::NO_TRANS,
                   Teuchos::ScalarTraits<scalar_t>::one(),
                   Teuchos::ScalarTraits<scalar_t>::zero());
}


AmrSmoother::Smoother
AmrSmoother::convertToEnumSmoother(const std::string& smoother) {
    std::map<std::string, Smoother> map;
    
    map["GS"]     = Smoother::GAUSS_SEIDEL;
    map["SGS"]    = Smoother::SGS;
    map["JACOBI"] = Smoother::JACOBI;
    
    auto sm = map.find(Util::toUpper(smoother));
    
    if ( sm == map.end() )
        throw OpalException("AmrMultiGrid::convertToEnumNorm_m()",
                            "No smoother '" + smoother + "'.");
    return sm->second;
}


void AmrSmoother::initParameter_m(const Smoother& smoother,
                                  lo_t nSweeps)
{
    if ( params_mp == Teuchos::null )
        params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    
    
    std::string type = "";
    scalar_t damping = 1.0;
    std::pair<bool, scalar_t> l1 = std::make_pair(true, 1.5);
    
    bool backward = false;
    std::pair<bool, scalar_t> fix = std::make_pair(true, 1.0e-5);
    bool check = true;
    
    switch ( smoother ) {
        
        case GAUSS_SEIDEL:
        {
            type = "Gauss-Seidel";
            backward = false;
            damping = 1.0;
            break;
        }
        case SGS:
        {
            type = "Symmetric Gauss-Seidel";
            damping = 1.0;
            break;
        }
        case JACOBI:
        {
            type = "Jacobi";
            damping = 6.0 / 7.0;
            break;
        }
        default:
            break;
    };
    
    params_mp->set("relaxation: type", type);
    params_mp->set("relaxation: sweeps", nSweeps);
    params_mp->set("relaxation: zero starting solution", false);
    params_mp->set("relaxation: damping factor", damping);
    params_mp->set("relaxation: use l1", l1.first);
    params_mp->set("relaxation: l1 eta", l1.second);
    params_mp->set("relaxation: backward mode", backward);
    params_mp->set("relaxation: fix tiny diagonal entries", fix.first);
    params_mp->set("relaxation: min diagonal value", fix.second);
    params_mp->set("relaxation: check diagonal entries", check);
}
