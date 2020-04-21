//
// Class AmrPoissonSolver
//   Base class of AMR poisson solvers.
//
// Copyright (c) 2016 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#ifndef AMR_POISSON_SOLVER_H_
#define AMR_POISSON_SOLVER_H_

#include "PoissonSolver.h"

#include <memory>

template<class AmrObject>
class AmrPoissonSolver : public PoissonSolver {
    
public:
    /*!
     * @param itsAmrObject_p holds information about grids and domain
     */
    AmrPoissonSolver(AmrObject* itsAmrObject_p)
        : itsAmrObject_mp(itsAmrObject_p), regrid_m(true) {}
    
    virtual ~AmrPoissonSolver() {}
    
    
    void computePotential(Field_t &rho, Vector_t hr) {
        throw OpalException("AmrPoissonSolver::computePotential(Field_t, Vector_t)",
                            "Not implemented.");
    }
    
    void computePotential(Field_t &rho, Vector_t hr, double zshift) {
        throw OpalException("AmrPoissonSolver::computePotential(Field_t, Vector_t, double)",
                            "Not implemented.");
    }
    
    void test(PartBunchBase<double, 3> *bunch) {
        throw OpalException("AmrPoissonSolver::test(PartBunchBase<double, 3>)", "Not implemented.");
    }
    
    void hasToRegrid() {
        regrid_m = true;
    }
    
protected:
    AmrObject* itsAmrObject_mp;
    
    /// is set to true by itsAmrObject_mp and reset to false by solver
    bool regrid_m;
};


#endif