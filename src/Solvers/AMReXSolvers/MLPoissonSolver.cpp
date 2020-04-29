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
#include "MLPoissonSolver.h"

#include "Utilities/OpalException.h"

#include "Amr/AmrDefs.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>

MLPoissonSolver::MLPoissonSolver(AmrBoxLib* itsAmrObject_p)
    : AmrPoissonSolver<AmrBoxLib>(itsAmrObject_p),
      reltol_m(1.0e-10),
      abstol_m(0.0)
{ }

void MLPoissonSolver::solve(AmrScalarFieldContainer_t& rho,
                            AmrScalarFieldContainer_t& phi,
                            AmrVectorFieldContainer_t& efield,
                            unsigned short baseLevel,
                            unsigned short finestLevel,
                            bool /*prevAsGuess*/)
{
    for (int i = 0; i <= finestLevel; ++i) {
        phi[i]->setVal(0.0, 1);
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield[i][j]->setVal(0.0, 1);
        }
    }
    
    this->mlmg_m(rho, phi, efield, baseLevel, finestLevel);
}


double MLPoissonSolver::getXRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(0);
}


double MLPoissonSolver::getXRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(0);
}


double MLPoissonSolver::getYRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(1);
}


double MLPoissonSolver::getYRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(1);
}


double MLPoissonSolver::getZRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(2);
}


double MLPoissonSolver::getZRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(2);
}


Inform &MLPoissonSolver::print(Inform &os) const {
    os << "* ************* A M R e X P o i s s o n S o l v e r ************************************ " << endl
       << "* relative tolerance " << reltol_m << '\n'
       << "* absolute tolerance " << abstol_m << '\n' << endl
       << "* ******************************************************************** " << endl;
    return os;
}


void MLPoissonSolver::mlmg_m(AmrScalarFieldContainer_t& rho,
                             AmrScalarFieldContainer_t& phi,
                             AmrVectorFieldContainer_t &efield,
                             int baseLevel,
                             int finestLevel)
{
    /* According to
     * amrex/Tutorials/LinearSolvers/ABecLaplacian_C/MyTest.H
     */
    amrex::LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMaxCoarseningLevel(10);
    
    const GeomContainer_t& geom = itsAmrObject_mp->Geom();
    
    int nlevels = finestLevel - baseLevel + 1;
    GeomContainer_t geom_p(nlevels);
    AmrGridContainer_t ba(nlevels);
    AmrProcMapContainer_t dm(nlevels);
    AmrFieldContainer_pt phi_p(nlevels);
    const_AmrFieldContainer_pt rho_p(nlevels);
    
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        geom_p[ilev] = geom[ilev + baseLevel];
        ba[ilev]     = rho[ilev + baseLevel]->boxArray();
        dm[ilev]     = rho[ilev + baseLevel]->DistributionMap();
        rho_p[ilev]  = rho[ilev + baseLevel].get();
        phi_p[ilev]  = phi[ilev + baseLevel].get();
    }
    
    amrex::MLPoisson mlpoisson(geom_p, ba, dm, info);
    
    mlpoisson.setMaxOrder(3);
    
    mlpoisson.setDomainBC({AMREX_D_DECL(amrex::LinOpBCType::Dirichlet,
                                        amrex::LinOpBCType::Dirichlet,
                                        amrex::LinOpBCType::Dirichlet)},
                          {AMREX_D_DECL(amrex::LinOpBCType::Dirichlet,
                                        amrex::LinOpBCType::Dirichlet,
                                        amrex::LinOpBCType::Dirichlet)});
    
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        mlpoisson.setLevelBC(ilev, phi[ilev].get());
    }
    
    amrex::MLMG mlmg(mlpoisson);
    mlmg.setMaxIter(200);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(0);
    mlmg.setCGVerbose(0);
    
    
    mlmg.solve(phi_p,
               rho_p,
               reltol_m, abstol_m);
    
    amrex::Vector<std::array<AmrField_t*, AMREX_SPACEDIM> > grad(nlevels);
    
    for (int i = 0; i < nlevels; ++i) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            grad[i][j] = efield[i][j].get();
        }
    }
    
    mlmg.getGradSolution(grad);
    
    for (int i = 0; i < nlevels; ++i) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield[i][j]->mult(-1.0, 0, 1, 1);
        }
    }
}
