//
// Class BottomSolver
//   Abstract base class for all base level solvers.
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
#ifndef BOTTOM_SOLVER_H
#define BOTTOM_SOLVER_H

#include "AmrMultiGridDefs.h"

template <class Matrix, class Vector, class Level>
class BottomSolver {
    
public:

    BottomSolver() : isInitialized_m(false) { };
    
    virtual ~BottomSolver() { };
    
    /*!
     * Solves
     * \f[
     *      Ax = b
     * \f]
     * @param x left-hand side
     * @param b right-hand side
     */
    virtual void solve(const Vector& x,
                       const Vector& b) = 0;
    
    /*!
     * Set the system matrix
     * @param A system matrix
     */
    virtual void setOperator(const Matrix& A,
                             Level* level_p = nullptr) = 0;
    
    
    /*!
     * @returns the number of required iterations
     */
    virtual std::size_t getNumIters() = 0;

    bool hasOperator() const;


protected:
    bool isInitialized_m;
};


template <class Matrix, class Vector, class Level>
bool BottomSolver<Matrix, Vector, Level>::hasOperator() const {
    return isInitialized_m;
}

#endif
