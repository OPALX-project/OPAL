//
// Class MLPoissonSolver
//   Interface to the C++ based AMR Poisson multigrid solver of AMReX.
//
// Copyright (c) 2016 - 2020, Matthias Frey, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved.
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
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef ML_POISSON_SOLVER_H_
#define ML_POISSON_SOLVER_H_

#include "Solvers/AmrPoissonSolver.h"
#include "Amr/AmrBoxLib.h"

class MLPoissonSolver : public AmrPoissonSolver< AmrBoxLib > {
    
private:
    typedef AmrBoxLib::AmrGeomContainer_t GeomContainer_t;
    typedef amrex::Vector<AmrBoxLib::AmrField_t*> AmrFieldContainer_pt;
    typedef amrex::Vector<const AmrBoxLib::AmrField_t*> const_AmrFieldContainer_pt;
    typedef AmrBoxLib::AmrGeometry_t AmrGeometry_t;
    typedef AmrBoxLib::AmrGrid_t AmrGrid_t;
    typedef AmrBoxLib::AmrProcMap_t           AmrProcMap_t;
    typedef AmrBoxLib::AmrScalarFieldContainer_t    AmrScalarFieldContainer_t;
    typedef AmrBoxLib::AmrVectorFieldContainer_t    AmrVectorFieldContainer_t;
    typedef AmrBoxLib::AmrGridContainer_t AmrGridContainer_t;
    typedef AmrBoxLib::AmrProcMapContainer_t AmrProcMapContainer_t;
    
public:
    
    /**
     * This solver only works with AmrBoxLib. In order the solver to work
     * the cells need to be of cubic shape. Otherwise the solver stops with
     * the error message that it did not converge (AMReX internal).
     * 
     * @param itsAmrObject_p has information about refinemen ratios, etc.
     */
    MLPoissonSolver(AmrBoxLib* itsAmrObject_p);
    
    /**
     * Multigrid solve based on AMReX FMultiGrid solver. The relative tolerance is
     * set to 1.0e-9 and the absolute tolerance to 0.0.
     * 
     * @param rho right-hand side charge density on grid [C / m]
     * @param phi electrostatic potential (unknown) [V]
     * @param efield electric field [V / m]
     * @param baseLevel for solve
     * @param finestLevel for solve
     * @param prevAsGuess use of previous solution as initial guess
     */
    void solve(AmrScalarFieldContainer_t& rho,
               AmrScalarFieldContainer_t& phi,
               AmrVectorFieldContainer_t& efield,
               unsigned short baseLevel,
               unsigned short finestLevel,
               bool prevAsGuess = true);
    
    double getXRangeMin(unsigned short level = 0);
    double getXRangeMax(unsigned short level = 0);
    double getYRangeMin(unsigned short level = 0);
    double getYRangeMax(unsigned short level = 0);
    double getZRangeMin(unsigned short level = 0);
    double getZRangeMax(unsigned short level = 0);
    
    
    /**
     * Print information abour tolerances.
     * @param os output stream where to write to
     */
    Inform &print(Inform &os) const;
    
    
private:
    
    /**
     * Does the actual solve.
     * 
     * @param rho charge density on grids [C / m]
     * @param phi electrostatic potential on grid [V]
     * @param baseLevel for solve
     * @param finestLevel for solve
     */
    void mlmg_m(AmrScalarFieldContainer_t& rho,
                AmrScalarFieldContainer_t& phi,
                AmrVectorFieldContainer_t &efield,
                int baseLevel,
                int finestLevel);
    
private:
    double reltol_m;                ///< Relative tolearance for solver
    double abstol_m;                ///< Absolute tolerance for solver
};


inline Inform &operator<<(Inform &os, const MLPoissonSolver &fs) {
    return fs.print(os);
}

#endif