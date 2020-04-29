//
// File AmrDefs
//   AMR namespace with typedefs of AMReX classes.
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
#ifndef AMR_DEFS_H
#define AMR_DEFS_H

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <memory>

/// Some AMR types used a lot
namespace amr {
    typedef amrex::MultiFab                                 AmrField_t;
    typedef std::array< std::unique_ptr<AmrField_t>,
                        AMREX_SPACEDIM
                       >                                    AmrVectorField_t;
    typedef amrex::DistributionMapping                      AmrProcMap_t;
    typedef amrex::Geometry                                 AmrGeometry_t;
    typedef amrex::BoxArray                                 AmrGrid_t;
    typedef amrex::Vector< std::unique_ptr<AmrField_t> >    AmrScalarFieldContainer_t;
    typedef amrex::Vector< AmrVectorField_t >               AmrVectorFieldContainer_t;
    typedef amrex::Vector< AmrGeometry_t >                  AmrGeomContainer_t;
    typedef amrex::Vector< AmrGrid_t >                      AmrGridContainer_t;
    typedef amrex::Vector< AmrProcMap_t >                   AmrProcMapContainer_t;
    typedef amrex::RealBox                                  AmrDomain_t;
    typedef amrex::Vector<int>                              AmrIntArray_t;
    typedef amrex::IntVect                                  AmrIntVect_t;
    typedef amrex::Vector< AmrIntVect_t >                   AmrIntVectContainer_t;
    typedef amrex::Box                                      AmrBox_t;
    typedef amrex::Real                                     AmrReal_t;
};

#endif
