#ifndef AMESOS2_SOLVER_H
#define AMESOS2_SOLVER_H

#include "BottomSolver.h"

#include <Amesos2.hpp>

#include <string>

extern Inform* gmsg;

/// Interface to Amesos2 solvers of the Trilinos package
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
