//
// Class Amesos2BottomSolver
//   Interface to Amesos2 solvers of the Trilinos package.
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
#ifndef AMESOS2_SOLVER_H
#define AMESOS2_SOLVER_H

#include "BottomSolver.h"

#include <Amesos2.hpp>

#include <string>

extern Inform* gmsg;

template <class Level>
class Amesos2BottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                                Teuchos::RCP<amr::multivector_t>,
                                                Level>
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::multivector_t mv_t;
    
    typedef Amesos2::Solver<matrix_t, mv_t> solver_t;
    
public:
    
    /*!
     * Instantiate
     * @param solvertype of Amesos2
     */
    Amesos2BottomSolver(std::string solvertype = "klu2");
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A,
                     Level* level_p = nullptr);
    
    std::size_t getNumIters();
    
private:
    
    std::string solvertype_m;           ///< kind of solver
    
    Teuchos::RCP<solver_t> solver_mp;   ///< solver instance
};

#include "Amesos2BottomSolver.hpp"

#endif
