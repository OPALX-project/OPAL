#ifndef BELOS_SOLVER_H
#define BELOS_SOLVER_H

#include "BottomSolver.h"

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
// #include <BelosSolverFactory_Tpetra.hpp>

#include "AmrPreconditioner.h"

#include <string>

/// Interface to Belos solvers of the Trilinos package
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
