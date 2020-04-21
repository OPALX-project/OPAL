//
// Class BelosBottomSolver
//   Interface to Belos solvers of the Trilinos package.
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
#ifndef BELOS_SOLVER_H
#define BELOS_SOLVER_H

#include "BottomSolver.h"

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

#include "AmrPreconditioner.h"

#include <string>

template <class Level>
class BelosBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                              Teuchos::RCP<amr::multivector_t>,
                                              Level>
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
    
    typedef Belos::SolverManager<scalar_t, mv_t, op_t> solver_t;
    typedef Belos::LinearProblem<scalar_t, mv_t, op_t> problem_t;
    
    typedef AmrPreconditioner<matrix_t, Level> prec_t;
    
public:
    
    /*!
     * @param solvertype to use
     * @param precond preconditioner of matrix
     */
    BelosBottomSolver(std::string solvertype = "Pseudoblock CG",
                      const std::shared_ptr<prec_t>& prec_p = nullptr);
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A,
                     Level* level_p = nullptr);
    
    std::size_t getNumIters();
    
private:
    /*!
     * Create a solver instance
     * @param solvertype to create
     */
    void initSolver_m(std::string solvertype);
    
private:
    Teuchos::RCP<problem_t> problem_mp;             ///< represents linear problem
    Teuchos::RCP<Teuchos::ParameterList> params_mp; ///< parameter list of solver
    Teuchos::RCP<solver_t>  solver_mp;              ///< solver instance
    std::shared_ptr<prec_t> prec_mp;                ///< preconditioner
    Teuchos::RCP<matrix_t> A_mp;                    ///< copy of matrix (has to be positive definite)
    
    scalar_t reltol_m;                              ///< relative tolerance
    
    /// allowed number of steps for iterative solvers
    int maxiter_m;
};

#include "BelosBottomSolver.hpp"

#endif
