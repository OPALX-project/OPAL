//
// Class BoxLibParticle
//   Particle class for AMReX. It works together with BoxLibLayout.
//   The class does the scatter and gather operations of attributes
//   to and from the grid. Ippl implements the same functionality in the
//   attribute class.
//
// Copyright (c) 2016 - 2020, Matthias Frey, Uldis Locans, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef BOXLIB_PARTICLE_HPP
#define BOXLIB_PARTICLE_HPP

#include "Utilities/OpalException.h"

#include <AMReX_BLFort.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_F.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

template<class PLayout>
BoxLibParticle<PLayout>::BoxLibParticle() : AmrParticleBase<PLayout>()
{
    AssignDensityTimer_m = IpplTimings::getTimer("AMR AssignDensity");
}


template<class PLayout>
BoxLibParticle<PLayout>::BoxLibParticle(PLayout *layout) : AmrParticleBase<PLayout>(layout)
{
    AssignDensityTimer_m = IpplTimings::getTimer("AMR AssignDensity");
}


template<class PLayout>
template<class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::scatter(ParticleAttrib<FT>& attrib, AmrScalarFieldContainer_t& f,
                                      ParticleAttrib<Vektor<PT, Dim> >& pp,
                                      int lbase, int lfine,
                                      const ParticleAttrib<int>& pbin, int bin)
{
    if ( lbase == lfine ) {
        this->scatter(attrib, *(f[lbase].get()), pp, pbin, bin, lbase);
        return;
    }
    
    const PLayout *layout_p = &this->getLayout();
    int nGrow = layout_p->refRatio(lbase)[0];
    
    AmrScalarFieldContainer_t tmp(lfine+1);
    for (int lev = lbase; lev <= lfine; ++lev) {
        
        f[lev]->setVal(0.0, f[lev]->nGrow());
        
        const AmrGrid_t& ba = f[lev]->boxArray();
        const AmrProcMap_t& dm = f[lev]->DistributionMap();
        tmp[lev].reset(new AmrField_t(ba, dm, 1, nGrow));
        tmp[lev]->setVal(0.0, nGrow);
    }
    
    this->AssignDensityFort(attrib, tmp, lbase, 1, lfine, pbin, bin);
    
    for (int lev = lbase; lev <= lfine; ++lev)
        AmrField_t::Copy(*f[lev], *tmp[lev], 0, 0, f[lev]->nComp(), f[lev]->nGrow());
}


template<class PLayout>
template <class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::scatter(ParticleAttrib<FT>& attrib, AmrField_t& f,
                                      ParticleAttrib<Vektor<PT, Dim> >& /*pp*/,
                                      const ParticleAttrib<int>& pbin, int bin,
                                      int level)
{
    const AmrGrid_t& ba      = f.boxArray();
    const AmrProcMap_t& dmap = f.DistributionMap();

    AmrField_t tmp(ba, dmap, f.nComp(), 1);
    tmp.setVal(0.0, 1);
    
    this->AssignCellDensitySingleLevelFort(attrib, tmp, level, pbin, bin);
    
    f.setVal(0.0, f.nGrow());
    
    AmrField_t::Copy(f, tmp, 0, 0, f.nComp(), f.nGrow());
}


template<class PLayout>
template<class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::gather(ParticleAttrib<FT>& attrib, AmrVectorFieldContainer_t& f,
                                     ParticleAttrib<Vektor<PT, Dim> >& /*pp*/,
                                     int lbase, int lfine)
{
    this->InterpolateFort(attrib, f, lbase, lfine);
}


// Function from BoxLib adjusted to work with Ippl BoxLibParticle class
// Scatter the particle attribute pa on the grid
template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::AssignDensityFort(ParticleAttrib<AType> &pa,
                                                AmrScalarFieldContainer_t& mf_to_be_filled, 
                                                int lev_min, int ncomp, int finest_level,
                                                const ParticleAttrib<int>& pbin, int bin) const
{
//     BL_PROFILE("AssignDensityFort()");
    IpplTimings::startTimer(AssignDensityTimer_m);
    const PLayout *layout_p = &this->getLayout();
    
    // not done in amrex
    int rho_index = 0;
    
    amrex::PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Vector<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;
    
    AmrScalarFieldContainer_t tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const AmrGrid_t& ba = mf_to_be_filled[lev]->boxArray();
        const AmrProcMap_t& dm = mf_to_be_filled[lev]->DistributionMap();
        tmp[lev].reset(new AmrField_t(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }
    
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        AssignCellDensitySingleLevelFort(pa, *mf_to_be_filled[lev], lev, pbin, bin, 1, 0);

        if (lev < finest_level) {
            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf_to_be_filled[lev],
                                          rho_index, rho_index, ncomp, 
                                          layout_p->Geom(lev), layout_p->Geom(lev+1),
                                          cphysbc, fphysbc,
                                          layout_p->refRatio(lev), &mapper, bcs);
        }

        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(*mf_to_be_filled[lev],
                                      *mf_to_be_filled[lev-1],
                                      rho_index, 1, layout_p->refRatio(lev-1),
                                      layout_p->Geom(lev-1), layout_p->Geom(lev));
        }
        
        mf_to_be_filled[lev]->plus(*tmp[lev], rho_index, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf_to_be_filled[lev+1], 
                             *mf_to_be_filled[lev], rho_index, ncomp, layout_p->refRatio(lev));
    }
    
    IpplTimings::stopTimer(AssignDensityTimer_m);
}


// This is the single-level version for cell-centered density
template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::AssignCellDensitySingleLevelFort(ParticleAttrib<AType> &pa,
                                                               AmrField_t& mf_to_be_filled,
                                                               int       level,
                                                               const     ParticleAttrib<int>& pbin,
                                                               int       bin,
                                                               int       ncomp,
                                                               int       /*particle_lvl_offset*/) const
{
//     BL_PROFILE("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignCellDensitySingleLevelFort()");
    
    const PLayout *layout_p = &this->getLayout();
    
    AmrField_t* mf_pointer;

    if (layout_p->OnSameGrids(level, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new AmrField_t(layout_p->ParticleBoxArray(level), 
                                  layout_p->ParticleDistributionMap(level),
                                  ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread 
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1)
        throw OpalException("BoxLibParticle::AssignCellDensitySingleLevelFort()",
                            "Must have at least one ghost cell when in AssignCellDensitySingleLevelFort");

    const AmrGeometry_t& gm          = layout_p->Geom(level);
    const AmrReal_t*     plo         = gm.ProbLo();
//     const AmrReal_t*     dx_particle = layout_p->Geom(level + particle_lvl_offset).CellSize();
    const AmrReal_t*     dx          = gm.CellSize();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
        throw OpalException("BoxLibParticle::AssignCellDensitySingleLevelFort()",
                            "Problem must be periodic in no or all directions");
    }
    
    for (amrex::MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    
    //loop trough particles and distribute values on the grid
    const ParticleLevelCounter_t& LocalNumPerLevel = this->getLocalNumPerLevel();
    size_t lBegin = LocalNumPerLevel.begin(level);
    size_t lEnd   = LocalNumPerLevel.end(level);
    
    AmrReal_t inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    
    for (size_t ip = lBegin; ip < lEnd; ++ip) {
        
        if ( bin > -1 && pbin[ip] != bin )
            continue;
        
        const int grid = this->Grid[ip];
        FArrayBox_t& fab = (*mf_pointer)[grid];
        
        // not callable:
        // begin amrex_deposit_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        AmrIntVect_t i1(i-1, j-1, k-1);
        AmrIntVect_t i2(i-1, j-1, k);
        AmrIntVect_t i3(i-1, j,   k-1);
        AmrIntVect_t i4(i-1, j,   k);
        AmrIntVect_t i5(i,   j-1, k-1);
        AmrIntVect_t i6(i,   j-1, k);
        AmrIntVect_t i7(i,   j,   k-1);
        AmrIntVect_t i8(i,   j,   k);
        
        fab(i1, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i2, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i3, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i4, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        fab(i5, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i6, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i7, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i8, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        // end of amrex_deposit_cic
    }
    
    mf_pointer->SumBoundary(gm.periodicity());

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const AmrReal_t vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }
}


template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::InterpolateFort(ParticleAttrib<AType> &pa,
                                              AmrVectorFieldContainer_t& mesh_data, 
                                              int lev_min, int lev_max)
{
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        if ( lev > 0 )
            InterpolateMultiLevelFort(pa, mesh_data, lev);
        else
            InterpolateSingleLevelFort(pa, mesh_data[lev], lev); 
    }
}


template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::InterpolateSingleLevelFort(ParticleAttrib<AType> &pa,
                                                         AmrVectorField_t& mesh_data, int lev)
{
    for (std::size_t i = 0; i < mesh_data.size(); ++i) {
        if (mesh_data[i]->nGrow() < 1)
            throw OpalException("BoxLibParticle::InterpolateSingleLevelFort()",
                                "Must have at least one ghost cell when in InterpolateSingleLevelFort");
    }
    
    PLayout *layout_p = &this->getLayout();
    
    const AmrGeometry_t& gm          = layout_p->Geom(lev);
    const AmrReal_t*     plo         = gm.ProbLo();
    const AmrReal_t*     dx          = gm.CellSize();
    
    //loop trough particles and distribute values on the grid
    const ParticleLevelCounter_t& LocalNumPerLevel = this->getLocalNumPerLevel();
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);    
    
    // make sure that boundaries are filled!
    for (std::size_t i = 0; i < mesh_data.size(); ++i) {
        mesh_data[i]->FillBoundary(gm.periodicity());
    }
    
    AmrReal_t inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    for (size_t ip = lBegin; ip < lEnd; ++ip) {

        const int grid = this->Grid[ip];
        FArrayBox_t& exfab = (*mesh_data[0])[grid];
        FArrayBox_t& eyfab = (*mesh_data[1])[grid];
        FArrayBox_t& ezfab = (*mesh_data[2])[grid];
        
        
        // not callable
        // begin amrex_interpolate_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), nComp, plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        AmrIntVect_t i1(i-1, j-1, k-1);
        AmrIntVect_t i2(i-1, j-1, k);
        AmrIntVect_t i3(i-1, j,   k-1);
        AmrIntVect_t i4(i-1, j,   k);
        AmrIntVect_t i5(i,   j-1, k-1);
        AmrIntVect_t i6(i,   j-1, k);
        AmrIntVect_t i7(i,   j,   k-1);
        AmrIntVect_t i8(i,   j,   k);
        
        pa[ip](0) = wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*exfab(i1) +
                    wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*exfab(i2) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*exfab(i3) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*exfab(i4) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*exfab(i5) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*exfab(i6) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*exfab(i7) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*exfab(i8);
        
        pa[ip](1) = wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*eyfab(i1) +
                    wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*eyfab(i2) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*eyfab(i3) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*eyfab(i4) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*eyfab(i5) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*eyfab(i6) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*eyfab(i7) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*eyfab(i8);
        
        pa[ip](2) = wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*ezfab(i1) +
                    wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*ezfab(i2) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*ezfab(i3) +
                    wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*ezfab(i4) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*ezfab(i5) +
                    wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*ezfab(i6) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*ezfab(i7) +
                    wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*ezfab(i8);
        // end amrex_interpolate_cic
    }
}


template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::InterpolateMultiLevelFort(ParticleAttrib<AType> &pa,
                                                        AmrVectorFieldContainer_t& mesh_data,
                                                        int lev)
{
    for (std::size_t i = 0; i < mesh_data[lev].size(); ++i) {
        if (mesh_data[lev][i]->nGrow() < 1)
            throw OpalException("BoxLibParticle::InterpolateMultiLevelFort()",
                                "Must have at least one ghost cell when in InterpolateMultiLevelFort");
    }
    
    PLayout *layout_p = &this->getLayout();
    
    const AmrGeometry_t& gm   = layout_p->Geom(lev);
    const AmrReal_t*     plo  = gm.ProbLo();
    const AmrReal_t*     fdx  = gm.CellSize();
    const AmrReal_t*     cdx  = layout_p->Geom(lev-1).CellSize();
    const AmrReal_t*     cplo = layout_p->Geom(lev-1).ProbLo();
    
    //loop trough particles and distribute values on the grid
    const ParticleLevelCounter_t& LocalNumPerLevel = this->getLocalNumPerLevel();
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);    
    
    // make sure that boundaries are filled!
    for (std::size_t i = 0; i < mesh_data[lev].size(); ++i) {
        mesh_data[lev][i]->FillBoundary(gm.periodicity());
    }
    
    AmrReal_t inv_fdx[3] = { 1.0 / fdx[0], 1.0 / fdx[1], 1.0 / fdx[2] };
    AmrReal_t inv_cdx[3] = { 1.0 / cdx[0], 1.0 / cdx[1], 1.0 / cdx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = { 0, 0, 0 };
    
    const AmrGrid_t& fba = mesh_data[lev][0]->boxArray();
    const AmrProcMap_t& fdmap = mesh_data[lev][0]->DistributionMap();
    AmrGrid_t cba = fba;
    cba.coarsen(AmrIntVect_t(2, 2, 2));

    AmrField_t cmesh_exdata(cba, fdmap,
                            mesh_data[lev][0]->nComp(),
                            mesh_data[lev][0]->nGrow());
    
    cmesh_exdata.setVal(0.0, 0, 1, mesh_data[lev][0]->nGrow());
    cmesh_exdata.copy(*mesh_data[lev-1][0], 0, 0, 1, 1, 1);
    cmesh_exdata.FillBoundary(gm.periodicity());
    
    AmrField_t cmesh_eydata(cba, fdmap,
                            mesh_data[lev][1]->nComp(),
                            mesh_data[lev][1]->nGrow());
    cmesh_eydata.setVal(0.0, 0, 1, mesh_data[lev][1]->nGrow());
    cmesh_eydata.copy(*mesh_data[lev-1][1], 0, 0, 1, 1, 1);
    cmesh_eydata.FillBoundary(gm.periodicity());
    
    AmrField_t cmesh_ezdata(cba, fdmap,
                            mesh_data[lev][2]->nComp(),
                            mesh_data[lev][2]->nGrow());
    cmesh_ezdata.setVal(0.0, 0, 1, mesh_data[lev][2]->nGrow());
    cmesh_ezdata.copy(*mesh_data[lev-1][2], 0, 0, 1, 1, 1);
    cmesh_ezdata.FillBoundary(gm.periodicity());

    for (size_t ip = lBegin; ip < lEnd; ++ip) {
        
        const int grid = this->Grid[ip];
        
        FArrayBox_t& exfab = (*(mesh_data[lev][0]))[grid];
        FArrayBox_t& eyfab = (*(mesh_data[lev][1]))[grid];
        FArrayBox_t& ezfab = (*(mesh_data[lev][2]))[grid];
        
        FArrayBox_t& cexfab = cmesh_exdata[grid];
        FArrayBox_t& ceyfab = cmesh_eydata[grid];
        FArrayBox_t& cezfab = cmesh_ezdata[grid];
        
        const typename PLayout::basefab_t& mfab = (*layout_p->getLevelMask(lev))[grid];
        
        // not callable
        // begin amrex_interpolate_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), nComp, plo, dx);
        for (int ii = 0; ii < 3; ++ii) {
            lxyz[ii] = ( this->R[ip](ii) - plo[ii] ) * inv_fdx[ii] + 0.5;
            ijk[ii] = lxyz[ii];
            wxyz_hi[ii] = lxyz[ii] - ijk[ii];
            wxyz_lo[ii] = 1.0 - wxyz_hi[ii];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        bool use_coarse = false;
        
        // AMReX: electrostatic_pic_3d.f90
        // use the coarse E near the level boundary
        if ( mfab(AmrIntVect_t(i-1, j-1, k-1)) == 1 ) {
            
            
            for (int ii = 0; ii < 3; ++ii) {
                lxyz[ii] = ( this->R[ip](ii) - cplo[ii] ) * inv_cdx[ii] + 0.5;
                ijk[ii] = lxyz[ii];
                wxyz_hi[ii] = lxyz[ii] - ijk[ii];
                wxyz_lo[ii] = 1.0 - wxyz_hi[ii];
            }
            use_coarse = true;
        }
        
        AmrIntVect_t i1(i-1, j-1, k-1);
        AmrIntVect_t i3(i-1, j,   k-1);
        AmrIntVect_t i5(i,   j-1, k-1);
        AmrIntVect_t i7(i,   j,   k-1);
        
        AmrIntVect_t i2(i-1, j-1, k);
        AmrIntVect_t i4(i-1, j,   k);
        AmrIntVect_t i6(i,   j-1, k);
        AmrIntVect_t i8(i,   j,   k);
        if ( use_coarse ) {
            pa[ip](0) = wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * cexfab(i1) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * cexfab(i3) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * cexfab(i5) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * cexfab(i7) +
                        wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * cexfab(i2) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * cexfab(i4) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * cexfab(i6) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * cexfab(i8);
                        
            pa[ip](1) = wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * ceyfab(i1) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * ceyfab(i3) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * ceyfab(i5) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * ceyfab(i7) +
                        wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * ceyfab(i2) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * ceyfab(i4) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * ceyfab(i6) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * ceyfab(i8);
                        
            pa[ip](2) = wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * cezfab(i1) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * cezfab(i3) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * cezfab(i5) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * cezfab(i7) +
                        wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * cezfab(i2) +
                        wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * cezfab(i4) +
                        wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * cezfab(i6) +
                        wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * cezfab(i8);
        } else {
            pa[ip](0) =  wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * exfab(i1) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * exfab(i3) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * exfab(i5) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * exfab(i7) +
                         wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * exfab(i2) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * exfab(i4) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * exfab(i6) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * exfab(i8);
            
            pa[ip](1) =  wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * eyfab(i1) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * eyfab(i3) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * eyfab(i5) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * eyfab(i7) +
                         wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * eyfab(i2) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * eyfab(i4) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * eyfab(i6) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * eyfab(i8);
            
            pa[ip](2) =  wxyz_lo[0] * wxyz_lo[1] * wxyz_lo[2] * ezfab(i1) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_lo[2] * ezfab(i3) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_lo[2] * ezfab(i5) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_lo[2] * ezfab(i7) +
                         wxyz_lo[0] * wxyz_lo[1] * wxyz_hi[2] * ezfab(i2) +
                         wxyz_lo[0] * wxyz_hi[1] * wxyz_hi[2] * ezfab(i4) +
                         wxyz_hi[0] * wxyz_lo[1] * wxyz_hi[2] * ezfab(i6) +
                         wxyz_hi[0] * wxyz_hi[1] * wxyz_hi[2] * ezfab(i8);
        }
        // end amrex_interpolate_cic
    }
}

#endif
