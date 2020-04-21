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
template <class Level>
Amesos2BottomSolver<Level>::Amesos2BottomSolver(std::string solvertype)
    : BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                Teuchos::RCP<amr::multivector_t>,
                                Level>()
    , solvertype_m(solvertype)
    , solver_mp(Teuchos::null)
{ }


template <class Level>
void Amesos2BottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                                       const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    solver_mp->solve(x.get(), b.get());
}


template <class Level>
void Amesos2BottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A,
                                             Level* level_p)
{
    try {
        solver_mp = Amesos2::create<matrix_t, mv_t>(solvertype_m, A);
    } catch(const std::invalid_argument& ex) {
        *gmsg << ex.what() << endl;
    }
    
    solver_mp->symbolicFactorization();
    solver_mp->numericFactorization();

    this->isInitialized_m = true;
}


template <class Level>
std::size_t Amesos2BottomSolver<Level>::getNumIters() {
    return 1;   // direct solvers do only one step
}
