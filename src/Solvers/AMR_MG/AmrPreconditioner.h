//
// Class AmrPreconditioner
//   Bottom solver preconditioners. Used with Belos bottom solvers.
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
#ifndef AMR_PRECONDITIONER_H
#define AMR_PRECONDITIONER_H

#include "AmrMultiGridDefs.h"

namespace amr {
    enum Preconditioner {
        NONE,
        ILUT,           ///< incomplete LU
        CHEBYSHEV,
        RILUK,          ///< ILU(k)
        SA,             ///< smoothed aggregation multigrid
        JACOBI,         ///< Jacobi point relaxation
        BLOCK_JACOBI,   ///< Jacobi block relaxation
        GS,             ///< Gauss-Seidel point relaxation
        BLOCK_GS        ///< Gauss-Seidel block relaxation
    };
}

template <class Matrix, class Level>
class AmrPreconditioner
{
public:
    typedef amr::operator_t operator_t;
    
public:
    
    /*!
     * Instantiate the preconditioner matrix
     * @param A matrix for which to create preconditioner
     * @param level_p bottom level if necessary to build preconditioner
     */
    virtual void create(const Teuchos::RCP<Matrix>& A, Level* level_p = nullptr) = 0;
    
    /*!
     * @returns the preconditioner
     */
    virtual Teuchos::RCP<operator_t> get() = 0;
};


#endif
