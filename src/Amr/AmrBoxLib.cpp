#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
#include "Structure/FieldSolver.h"
#include "Solvers/PoissonSolver.h"

#include "Amr/AmrYtWriter.h"
    
#include <AMReX_MultiFabUtil.H>

#include <AMReX_ParmParse.H> // used in initialize function
#include <AMReX_BCUtil.H>

#include <map>

extern Inform* gmsg;


AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     int maxLevel,
                     AmrPartBunch* bunch_p)
    : AmrObject()
    , amrex::AmrMesh(&domain, maxLevel, nGridPts, 0 /* cartesian */)
    , bunch_mp(bunch_p)
    , layout_mp(static_cast<AmrLayout_t*>(&bunch_p->getLayout()))
    , rho_m(maxLevel + 1)
    , phi_m(maxLevel + 1)
    , efield_m(maxLevel + 1)
    , meshScaling_m(Vector_t(1.0, 1.0, 1.0))
    , isFirstTagging_m(maxLevel + 1, true)
    , isPoissonSolved_m(false)
{
    /*
     * The layout needs to know how many levels we can make.
     */
    layout_mp->resize(maxLevel);
    
    initBaseLevel_m(nGridPts);
    
    // set mesh spacing of bunch
    updateMesh();
}


std::unique_ptr<AmrBoxLib> AmrBoxLib::create(const AmrInfo& info,
                                             AmrPartBunch* bunch_p)
{
    /* The bunch is initialized first with a Geometry,
     * BoxArray and DistributionMapping on
     * the base level (level = 0). Thus, we take the domain specified there in
     * order to create the Amr object.
     */
    AmrLayout_t* layout_p = static_cast<AmrLayout_t*>(&bunch_p->getLayout());
    AmrDomain_t domain = layout_p->Geom(0).ProbDomain();
    
    AmrIntArray_t nGridPts = {
        info.grid[0],
        info.grid[1],
        info.grid[2]
    };
    
    int maxlevel = info.maxlevel;
    
    /*
     * further attributes are given by the BoxLib's ParmParse class.
     */
    initParmParse_m(info, layout_p);
    
    return std::unique_ptr<AmrBoxLib>(new AmrBoxLib(domain,
                                                    nGridPts,
                                                    maxlevel,
                                                    bunch_p
                                                   )
                                     );
}


void AmrBoxLib::regrid(double time) {
    IpplTimings::startTimer(this->amrRegridTimer_m);

    *gmsg << "* Start regriding:" << endl
          << "*     Old finest level: "
          << finest_level << endl;

    /* ATTENTION: The bunch might be updated during
     * the regrid process!
     * We regrid from base level 0 up to the finest level.
     */
    int old_finest = finest_level;
    
    int lev_top = std::min(finest_level, max_level - 1);
    for (int i = 0; i <= lev_top; ++i) {
        this->doRegrid_m(i, time);
        lev_top = std::min(finest_level, max_level - 1);
    }

    this->postRegrid_m(old_finest);

    *gmsg << "*     New finest level: "
          << finest_level << endl
          << "* Finished regriding" << endl;

    IpplTimings::stopTimer(this->amrRegridTimer_m);
}


void AmrBoxLib::getGridStatistics(std::map<int, long>& gridPtsPerCore,
                                  std::vector<int>& gridsPerLevel) const
{
    typedef std::vector<int> container_t;
    
    gridPtsPerCore.clear();
    gridsPerLevel.clear();
    
    gridsPerLevel.resize(max_level + 1);
    
    for (int lev = 0; lev <= finest_level; ++lev) {
        /* container index: box
         * container value: cores that owns box
         */
        const container_t& pmap = this->dmap[lev].ProcessorMap();
        const AmrGrid_t& ba = this->grids[lev];
        
        gridsPerLevel[lev] = pmap.size();

        // iterate over all boxes
        for (unsigned int i = 0; i < ba.size(); ++i) {
            gridPtsPerCore[pmap[i]] += ba[i].numPts();
        }
    }
}


void AmrBoxLib::initFineLevels() {
    if ( bunch_mp->getNumBunch() > 1 && !refined_m ) {
        *gmsg << "* Initialization of all levels" << endl;
        
        AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
        
        /* we do an explicit domain mapping of the particles and then
         * forbid it during the regrid process, this way it's only
         * executed ones --> saves computation
         */
        bool isForbidTransform = amrpbase_p->isForbidTransform();
        
        if ( !isForbidTransform ) {
            amrpbase_p->domainMapping();
            amrpbase_p->setForbidTransform(true);
        }

        if ( max_level > 0 ) {
            this->regrid(bunch_mp->getT() * 1.0e9 /*time [ns] */);
        }

        if ( !isForbidTransform ) {
            amrpbase_p->setForbidTransform(false);
            // map particles back
            amrpbase_p->domainMapping(true);
        }
        
        *gmsg << "* Initialization done." << endl;
        
        refined_m = true;
    }
}


AmrBoxLib::VectorPair_t AmrBoxLib::getEExtrema() {
    Vector_t maxE = Vector_t(0, 0, 0), minE = Vector_t(0, 0, 0);
    /* check for maximum / minimum values over all levels
     * and all components, i.e. Ex, Ey, Ez
     */
    for (unsigned int lev = 0; lev < efield_m.size(); ++lev) {
        for (std::size_t i = 0; i < bunch_mp->Dimension; ++i) {
            /* calls:
             *  - max(comp, nghost (to search), local)
             *  - min(cmop, nghost, local)
             */
            double max = efield_m[lev][i]->max(0, false);
            maxE[i] = (maxE[i] < max) ? max : maxE[i];
            
            double min = efield_m[lev][i]->min(0, false);
            minE[i] = (minE[i] > min) ? min : minE[i];
        }
    }
    
    return VectorPair_t(maxE, minE);
}


double AmrBoxLib::getRho(int x, int y, int z) {
    //TODO
    throw OpalException("AmrBoxLib::getRho(x, y, z)", "Not yet Implemented.");
    return 0.0;
}


void AmrBoxLib::computeSelfFields() {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields(int bin)", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields_cycl(double gamma) {   
    /*
     * The potential is not scaled according to domain modification.
     * 
     */
    if ( !bunch_mp->hasFieldSolver() )
        return;
    
    /*
     * scatter charges onto grid
     */
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    
    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
    bunch_mp->updateLorentzFactor(gamma);
    
    // map on Amr domain + Lorentz transform
    double scalefactor = amrpbase_p->domainMapping();
    
    amrpbase_p->setForbidTransform(true);
    amrpbase_p->update();
    amrpbase_p->setForbidTransform(false);
    
    /// from charge (C) to charge density (C/m^3).
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level, bunch_mp->Bin);
    
    int baseLevel = 0;
    int nLevel = finest_level + 1;
    double invGamma = 1.0 / gamma;
    
    
    // mesh scaling for solver
    meshScaling_m = Vector_t(1.0, 1.0, 1.0);
    
    // charge density is in rho_m
    // calculate Possion equation (with coefficient: -1/(eps))
    for (int i = 0; i <= finest_level; ++i) {
        if ( this->rho_m[i]->contains_nan(false) )
            throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                "NANs at level " + std::to_string(i) + ".");
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    }
    
    // find maximum and normalize each level (faster convergence)
    double l0norm = 1.0;
    for (int i = 0; i <= finest_level; ++i)
        l0norm = std::max(l0norm, this->rho_m[i]->norm0(0));
    
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(1.0 / l0norm, 0, 1);
    }

    PoissonSolver *solver = bunch_mp->getFieldSolver();

    IpplTimings::startTimer(this->amrSolveTimer_m);
    solver->solve(rho_m, phi_m, efield_m, baseLevel, finest_level);
    IpplTimings::stopTimer(this->amrSolveTimer_m);
    
    // make sure ghost cells are filled
    for (int i = 0; i <= finest_level; ++i) {
        phi_m[i]->FillBoundary(geom[i].periodicity());
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield_m[i][j]->FillBoundary(geom[i].periodicity());
        }
    }
    
    this->fillPhysbc_m(*(this->phi_m[0]), 0);
    
    /* apply scale of electric-field in order to undo the transformation
     * + undo normalization
     */
    for (int i = 0; i <= finestLevel(); ++i) {
        this->phi_m[i]->mult(scalefactor * l0norm, 0, 1);
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            this->efield_m[i][j]->mult(scalefactor * scalefactor * l0norm, 0, 1);
        }
    }
    
    for (int i = 0; i <= finest_level; ++i) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            if ( this->efield_m[i][j]->contains_nan(false) )
                throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                    "Ef: NANs at level " + std::to_string(i) + ".");
        }
    }
    
    
    amrpbase_p->gather(bunch_mp->Ef, this->efield_m, bunch_mp->R, 0, finest_level);
    
    // undo domain change + undo Lorentz transform
    amrpbase_p->domainMapping(true);
    
    /// Back Lorentz transformation
    bunch_mp->Ef *= Vector_t(gamma,
                             1.0,
                             gamma);
    
    /// calculate coefficient
    // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
    // but because we already transformed E_trans into the moving frame we have to
    // add 1/gamma so we are using the E_trans from the rest frame -DW
    double betaC = std::sqrt(gamma * gamma - 1.0) * invGamma / Physics::c;
    
    /// calculate B field from E field
    bunch_mp->Bf(0) =  betaC * bunch_mp->Ef(2);
    bunch_mp->Bf(2) = -betaC * bunch_mp->Ef(0);
    
    /*
     * dumping only
     */
    
    if ( !(bunch_mp->getLocalTrackStep()  % Options::amrYtDumpFreq) ) {
        AmrYtWriter ytWriter(bunch_mp->getLocalTrackStep());
    
        AmrIntArray_t rr(finest_level);
        for (int i = 0; i < finest_level; ++i)
            rr[i] = this->MaxRefRatio(i);
        
        double time = bunch_mp->getT(); // in seconds
        
        // we need to undo coefficient when writing charge density
        for (int i = 0; i <= finest_level; ++i)
            this->rho_m[i]->mult(- Physics::epsilon_0 * l0norm, 0, 1);
        
        
        ytWriter.writeFields(rho_m, phi_m, efield_m, rr, this->geom,
                             nLevel, time, scalefactor);
        
        ytWriter.writeBunch(bunch_mp, time, gamma);
    }
}


void AmrBoxLib::computeSelfFields_cycl(int bin) {

    if ( !bunch_mp->hasFieldSolver() )
        return;
    
    /// get gamma of this bin
    double gamma = bunch_mp->getBinGamma(bin);
    
    /* relativistic factor is always gamma >= 1
     * --> if no particle --> gamma = 0 --> leave computation
     */
    if ( gamma < 1.0 )
        return;

    /*
     * scatter charges onto grid
     */
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    
    // Lorentz transformation: apply Lorentz factor of that particle bin
    bunch_mp->updateLorentzFactor(bin);
    
    // map on Amr domain + Lorentz transform 
    double scalefactor = amrpbase_p->domainMapping();
    
    amrpbase_p->setForbidTransform(true);
    
    amrpbase_p->update();

    if ( !(bunch_mp->getLocalTrackStep() % Options::amrRegridFreq) ) {
        this->regrid(bunch_mp->getT() * 1.0e9 /*time [ns] */);
    }

    amrpbase_p->setForbidTransform(false);
    
    /// scatter particles charge onto grid.
    /// from charge (C) to charge density (C/m^3).
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level, bunch_mp->Bin, bin);
    
    /// Lorentz transformation
    // In particle rest frame, the longitudinal length (y for cyclotron) enlarged
    // calculate Possion equation (with coefficient: -1/(eps))
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(-1.0 / (Physics::epsilon_0), 0 /*comp*/, 1 /*ncomp*/);
        
        if ( this->rho_m[i]->contains_nan(false) )
            throw OpalException("AmrBoxLib::computeSelfFields_cycl(int bin) ",
                                "NANs at level " + std::to_string(i) + ".");
    }
    
    // find maximum and normalize each level (faster convergence)
    double l0norm = 1.0;
    for (int i = 0; i <= finest_level; ++i)
        l0norm = std::max(l0norm, this->rho_m[i]->norm0(0));
    
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    // mesh scaling for solver
    meshScaling_m = Vector_t(1.0, 1.0, 1.0);
    
    PoissonSolver *solver = bunch_mp->getFieldSolver();

    IpplTimings::startTimer(this->amrSolveTimer_m);
    
    for (int i = 0; i <= finest_level; ++i) {
        phi_m[i]->setVal(0.0, 0, phi_m[i]->nComp(), phi_m[i]->nGrow());
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield_m[i][j]->setVal(0.0, 0, efield_m[i][j]->nComp(), efield_m[i][j]->nGrow());
        }
    }
    
    // in case of binning we reset phi every time
    solver->solve(rho_m, phi_m, efield_m, 0, finest_level, false);
    
    IpplTimings::stopTimer(this->amrSolveTimer_m);
    
    // make sure ghost cells are filled
    for (int i = 0; i <= finest_level; ++i) {
        phi_m[i]->FillBoundary(geom[i].periodicity());
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            efield_m[i][j]->FillBoundary(geom[i].periodicity());
        }
    }
    
    this->fillPhysbc_m(*(this->phi_m[0]), 0);
    
    
    /* apply scale of electric-field in order to undo the transformation
     * + undo normalization
     */
    for (int i = 0; i <= finestLevel(); ++i) {
        this->phi_m[i]->mult(scalefactor * l0norm, 0, 1);
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            this->efield_m[i][j]->mult(scalefactor * scalefactor * l0norm, 0, 1);
        }
    }

    for (int i = 0; i <= finest_level; ++i) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            if ( this->efield_m[i][j]->contains_nan(false) )
                throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                    "Ef: NANs at level " + std::to_string(i) + ".");
        }
    }
    
    amrpbase_p->gather(bunch_mp->Eftmp, this->efield_m, bunch_mp->R, 0, finest_level);
    
    // undo domain change + undo Lorentz transform
    amrpbase_p->domainMapping(true);
    
    /// Back Lorentz transformation
    bunch_mp->Eftmp *= Vector_t(gamma, 1.0, gamma);

    /// Calculate coefficient
    double betaC = std::sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

    /// Calculate B_bin field from E_bin field accumulate B and E field
    bunch_mp->Bf(0) += betaC * bunch_mp->Eftmp(2);
    bunch_mp->Bf(2) -= betaC * bunch_mp->Eftmp(0);

    bunch_mp->Ef += bunch_mp->Eftmp;

    /*
     * dumping only
     */
    if ( !(bunch_mp->getLocalTrackStep()  % Options::amrYtDumpFreq) ) {
        AmrYtWriter ytWriter(bunch_mp->getLocalTrackStep(), bin);
        
        int nLevel = finest_level + 1;
        
        AmrIntArray_t rr(finest_level);
        for (int i = 0; i < finest_level; ++i)
            rr[i] = this->MaxRefRatio(i);
        
        double time = bunch_mp->getT(); // in seconds
        
        // we need to undo coefficient when writing charge density
        for (int i = 0; i <= finest_level; ++i)
            this->rho_m[i]->mult(- Physics::epsilon_0 * l0norm, 0, 1);
        
        
        ytWriter.writeFields(rho_m, phi_m, efield_m, rr, this->geom,
                             nLevel, time, scalefactor);
        
        ytWriter.writeBunch(bunch_mp, time, gamma);
    }

    isPoissonSolved_m = true;
}


void AmrBoxLib::updateMesh() {
    //FIXME What about resizing mesh, i.e. geometry?    
    const AmrReal_t* tmp = this->geom[0].CellSize();
    
    Vector_t hr;
    for (int i = 0; i < 3; ++i)
        hr[i] = tmp[i];
    
    bunch_mp->setBaseLevelMeshSpacing(hr);
}


const Vector_t& AmrBoxLib::getMeshScaling() const {
    return meshScaling_m;
}


Vektor<int, 3> AmrBoxLib::getBaseLevelGridPoints() const {
    const Box_t& bx = this->geom[0].Domain();
    
    const AmrIntVect_t& low = bx.smallEnd();
    const AmrIntVect_t& high = bx.bigEnd();
    
    return Vektor<int, 3>(high[0] - low[0] + 1,
                          high[1] - low[1] + 1,
                          high[2] - low[2] + 1);
}


inline const int& AmrBoxLib::maxLevel() const {
    return this->max_level;
}


inline const int& AmrBoxLib::finestLevel() const {
    return this->finest_level;
}


inline double AmrBoxLib::getT() const {
    return bunch_mp->getT();
}


void AmrBoxLib::redistributeGrids(int how) {
//
//    // copied + modified version of AMReX_Amr.cpp
//    AmrProcMap_t::InitProximityMap();
////     AmrProcMap_t::Initialize();
//
//    AmrGridContainer_t allBoxes(finest_level + 1);
//    for(unsigned int ilev = 0; ilev < allBoxes.size(); ++ilev) {
//        allBoxes[ilev] = boxArray(ilev);
//    }
//    amrex::Vector<AmrIntArray_t> mLDM;
//        
//    switch ( how ) {
//        case RANK_ZERO:
//            mLDM = AmrProcMap_t::MultiLevelMapRandom(this->ref_ratio,
//                                                     allBoxes,
//                                                     maxGridSize(0),
//                                                     0/*maxRank*/,
//                                                     0/*minRank*/);
//            break;
//        case PFC:
//            mLDM = AmrProcMap_t::MultiLevelMapPFC(this->ref_ratio,
//                                                  allBoxes,
//                                                  maxGridSize(0));
//            break;
//        case RANDOM:
//            mLDM = AmrProcMap_t::MultiLevelMapRandom(this->ref_ratio,
//                                                     allBoxes,
//                                                     maxGridSize(0));
//            break;
//        case KNAPSACK:
//            mLDM = AmrProcMap_t::MultiLevelMapKnapSack(this->ref_ratio,
//                                                       allBoxes,
//                                                       maxGridSize(0));
//            break;
//        default:
//            *gmsg << "We didn't redistribute the grids." << endl;
//            return;
//    }
//        
//    for(unsigned int iMap = 0; iMap < mLDM.size(); ++iMap) {
//        AmrField_t::MoveAllFabs(mLDM[iMap]);
//    }
//    
//    /*
//     * particles need to know the BoxArray
//     * and DistributionMapping
//     */
//    for(unsigned int ilev = 0; ilev < allBoxes.size(); ++ilev) {
//        layout_mp->SetParticleBoxArray(ilev, this->grids[ilev]);
//        layout_mp->SetParticleDistributionMap(ilev, this->dmap[ilev]);
//    }
//    
//    for(unsigned int iMap = 0; iMap < mLDM.size(); ++iMap) {
//        rho_m[iMap]->MoveFabs(mLDM[iMap]);
//        phi_m[iMap]->MoveFabs(mLDM[iMap]);
//        efield_m[iMap]->MoveFabs(mLDM[iMap]);
//    }
}


void AmrBoxLib::RemakeLevel (int lev, AmrReal_t time,
                             const AmrGrid_t& new_grids,
                             const AmrProcMap_t& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    //                                                      #comp  #ghosts cells
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     0));
    phi_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        efield_m[lev][j].reset(new AmrField_t(new_grids, new_dmap, 1,     1));
    }
    
    // including nghost = 1
    rho_m[lev]->setVal(0.0, 0);
    phi_m[lev]->setVal(0.0, 1);
    
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        efield_m[lev][j]->setVal(0.0, 1);
    }
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
     layout_mp->SetParticleBoxArray(lev, new_grids);
     layout_mp->SetParticleDistributionMap(lev, new_dmap);
     
     layout_mp->buildLevelMask(lev, this->nErrorBuf(lev));
}


void AmrBoxLib::MakeNewLevel (int lev, AmrReal_t time,
                              const AmrGrid_t& new_grids,
                              const AmrProcMap_t& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    //                                                      #comp  #ghosts cells
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     0));
    phi_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        efield_m[lev][j].reset(new AmrField_t(new_grids, new_dmap, 1,     1));
    }
    
    // including nghost = 1
    rho_m[lev]->setVal(0.0, 0);
    phi_m[lev]->setVal(0.0, 1);
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        efield_m[lev][j]->setVal(0.0, 1);
    }
    
    
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
    layout_mp->SetParticleBoxArray(lev, new_grids);
    layout_mp->SetParticleDistributionMap(lev, new_dmap);
    
    layout_mp->buildLevelMask(lev, this->nErrorBuf(lev));
}


void AmrBoxLib::ClearLevel(int lev) {
    rho_m[lev].reset(nullptr);
    phi_m[lev].reset(nullptr);
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        efield_m[lev][j].reset(nullptr);
    }
    ClearBoxArray(lev);
    ClearDistributionMap(lev);

    this->isFirstTagging_m[lev] = true;
}


void AmrBoxLib::ErrorEst(int lev, TagBoxArray_t& tags,
                         AmrReal_t time, int ngrow)
{
    *gmsg << level2 << "*         Start tagging of level " << lev << endl;
    
    switch ( tagging_m ) {
        case CHARGE_DENSITY:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
        case POTENTIAL:
            tagForPotentialStrength_m(lev, tags, time, ngrow);
            break;
        case EFIELD:
            tagForEfield_m(lev, tags, time, ngrow);
            break;
        case MOMENTA:
            tagForMomenta_m(lev, tags, time, ngrow);
            break;
        case MAX_NUM_PARTICLES:
            tagForMaxNumParticles_m(lev, tags, time, ngrow);
            break;
        case MIN_NUM_PARTICLES:
            tagForMinNumParticles_m(lev, tags, time, ngrow);
            break;
        default:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
    }
    
    *gmsg << level2 << "*         Finished tagging of level " << lev << endl;
}


void AmrBoxLib::MakeNewLevelFromScratch(int lev, AmrReal_t time,
                                  const AmrGrid_t& ba,
                                  const AmrProcMap_t& dm)
{
    throw OpalException("AmrBoxLib::MakeNewLevelFromScratch()", "Shouldn't be called.");
}


void AmrBoxLib::MakeNewLevelFromCoarse (int lev, AmrReal_t time,
                                        const AmrGrid_t& ba,
                                        const AmrProcMap_t& dm)
{
    throw OpalException("AmrBoxLib::MakeNewLevelFromCoarse()", "Shouldn't be called.");
}


void AmrBoxLib::doRegrid_m(int lbase, double time) {
    int new_finest = 0;
    AmrGridContainer_t new_grids(finest_level+2);
    
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);
    
    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
        if (lev <= finest_level) // an old level
        {
            if (new_grids[lev] != grids[lev]) // otherwise nothing
            {
                AmrProcMap_t new_dmap(new_grids[lev]);
                RemakeLevel(lev, time, new_grids[lev], new_dmap);
            }
        }
        else  // a new level
        {
            AmrProcMap_t new_dmap(new_grids[lev]);
            MakeNewLevel(lev, time, new_grids[lev], new_dmap);
        }
        
        layout_mp->setFinestLevel(new_finest);
    }

    // now we are safe to delete deprecated levels
    for (int lev = new_finest+1; lev <= finest_level; ++lev) {
        ClearLevel(lev);
    }
    
    finest_level = new_finest;
}


void AmrBoxLib::postRegrid_m(int old_finest) {
    /* ATTENTION: In this call the bunch has to be updated
     * since each particle needs to be assigned to its latest
     * grid and sent to the corresponding MPI-process.
     * 
     * We need to update the bunch before we clear
     * levels, otherwise the particle level counter
     * is invalidated.
     */

    // redistribute particles and assign correct grid and level
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->update(0, finest_level, true);

    // check if no particles left on higher levels
    const auto& LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();

    for (int lev = finest_level+1; lev <= old_finest; ++lev) {
        if ( LocalNumPerLevel.getLocalNumAtLevel(lev) != 0 ) {
            throw OpalException("AmrBoxLib::postRegrid_m()",
                                "Still particles on level " + std::to_string(lev) + "!");
        }
        layout_mp->ClearParticleBoxArray(lev);
        layout_mp->ClearParticleDistributionMap(lev);
        layout_mp->clearLevelMask(lev);
    }

    if ( bunch_mp->getLocalNum() != LocalNumPerLevel.getLocalNumUpToLevel(finest_level) ) {
        std::string localnum = std::to_string(bunch_mp->getLocalNum());
        std::string levelnum = std::to_string(LocalNumPerLevel.getLocalNumUpToLevel(finest_level));
        throw OpalException("AmrBoxLib::postRegrid_m()",
                            "Number of particles do not agree: " + localnum + " != " + levelnum);
    }

    PoissonSolver *solver = bunch_mp->getFieldSolver();
    solver->hasToRegrid();
}


void AmrBoxLib::tagForChargeDensity_m(int lev, TagBoxArray_t& tags,
                                      AmrReal_t time, int ngrow)
{
    // we need to update first
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->update(0, lev, true);
    
    for (int i = lev; i <= finest_level; ++i)
        rho_m[i]->setVal(0.0, rho_m[i]->nGrow());
    
    // the new scatter function averages the value also down to the coarsest level
    amrpbase_p->scatter(bunch_mp->Q, rho_m,
                        bunch_mp->R, 0, lev, bunch_mp->Bin);
    
    const double& scalefactor = amrpbase_p->getScalingFactor();
    
    for (int i = lev; i <= finest_level; ++i)
        rho_m[i]->mult(scalefactor, 0, 1);
    
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(*rho_m[lev], true); mfi.isValid(); ++mfi) {
            const Box_t&  tilebx  = mfi.tilebox();
            
            TagBox_t&     tagfab  = tags[mfi];
            FArrayBox_t&  fab = (*rho_m[lev])[mfi];
            
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();
            
            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        
                        amrex::IntVect iv(D_DECL(i,j,k));
                        
                        if ( std::abs( fab(iv) ) >= chargedensity_m )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::tagForPotentialStrength_m(int lev, TagBoxArray_t& tags,
                                          AmrReal_t time, int ngrow)
{
    /* Tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    if ( !isPoissonSolved_m || !phi_m[lev]->ok() || this->isFirstTagging_m[lev] ) {
        *gmsg << level2 << "* Level " << lev << ": We need to perform "
              << "charge tagging if a new level is created." << endl;
        this->tagForChargeDensity_m(lev, tags, time, ngrow);
        this->isFirstTagging_m[lev] = false;
        return;
    }
    
    // tag cells for refinement
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;

    AmrReal_t threshold = phi_m[lev]->norm0(0, 0 /*nghost*/, false /*local*/);
    
    threshold *= scaling_m;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(*phi_m[lev], true); mfi.isValid(); ++mfi) {
            
            const Box_t&  tilebx  = mfi.tilebox();
            TagBox_t&     tagfab  = tags[mfi];
            FArrayBox_t&  fab     = (*phi_m[lev])[mfi];
            
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();
            
            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        
                        amrex::IntVect iv(D_DECL(i,j,k));
                        
                        if ( std::abs( fab(iv) ) >= threshold )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray_t& tags,
                               AmrReal_t time, int ngrow)
{
    /* Tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    if ( !isPoissonSolved_m || !efield_m[lev][0]->ok() || this->isFirstTagging_m[lev] ) {
        *gmsg << level2 << "* Level " << lev << ": We need to perform "
              << "charge tagging if a new level is created." << endl;
        this->tagForChargeDensity_m(lev, tags, time, ngrow);
        this->isFirstTagging_m[lev] = false;
        return;
    }
    
    // tag cells for refinement
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;
    
    // obtain maximum absolute value
    amrex::Vector<AmrReal_t> threshold(AMREX_SPACEDIM);
    
    for (int i = 0; i < 3; ++i) {
        threshold[i] = efield_m[lev][i]->norm0(0,
                                            0 /*nghosts*/,
                                            false /*local*/);
    
        threshold[i] *= scaling_m;
    }
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(this->grids[lev], this->dmap[lev], true); mfi.isValid(); ++mfi) {
            
            const Box_t&  tilebx  = mfi.tilebox();
            TagBox_t&     tagfab  = tags[mfi];
            FArrayBox_t&  xfab     = (*efield_m[lev][0])[mfi];
            FArrayBox_t&  yfab     = (*efield_m[lev][1])[mfi];
            FArrayBox_t&  zfab     = (*efield_m[lev][2])[mfi];

            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        AmrIntVect_t iv(i,j,k);
                        if (std::abs(xfab(iv)) >= threshold[0])
                            tagfab(iv) = tagval;
                        else if (std::abs(yfab(iv)) >= threshold[1])
                            tagfab(iv) = tagval;
                        else if (std::abs(zfab(iv)) >= threshold[2])
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::tagForMomenta_m(int lev, TagBoxArray_t& tags,
                                AmrReal_t time, int ngrow)
{
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // we need to update first
    amrpbase_p->update(0, lev, true);
    const auto& LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();
    
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);
    
    Vector_t pmax = Vector_t(0.0, 0.0, 0.0);
    for (size_t i = lBegin; i < lEnd; ++i) {
        const Vector_t& tmp = bunch_mp->P[i];
        pmax = ( dot(tmp, tmp) > dot(pmax, pmax) ) ? tmp : pmax;
    }
    
    double momentum2 = scaling_m * dot(pmax, pmax);
    
    std::map<AmrIntVect_t, bool> cells;
    for (size_t i = lBegin; i < lEnd; ++i) {
        if ( dot(bunch_mp->P[i], bunch_mp->P[i]) >= momentum2 ) {
            AmrIntVect_t iv = layout_mp->Index(bunch_mp->R[i], lev);
            cells[iv] = true;
        }
    }
        
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(*rho_m[lev], true); mfi.isValid(); ++mfi) {
            
            const Box_t&  tilebx  = mfi.tilebox();
            TagBox_t&     tagfab  = tags[mfi];
            
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        AmrIntVect_t iv(i, j, k);
                        if ( cells.find(iv) != cells.end() )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::tagForMaxNumParticles_m(int lev, TagBoxArray_t& tags,
                                        AmrReal_t time, int ngrow)
{
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // we need to update first
    amrpbase_p->update(0, lev, true);
    const auto& LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();
    
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);
    
    // count particles per cell
    std::map<AmrIntVect_t, size_t> cells;
    for (size_t i = lBegin; i < lEnd; ++i) {
        AmrIntVect_t iv = layout_mp->Index(bunch_mp->R[i], lev);
        ++cells[iv];
    }
        
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(*rho_m[lev], true); mfi.isValid(); ++mfi) {
            
            const Box_t&  tilebx  = mfi.tilebox();
            TagBox_t&     tagfab  = tags[mfi];

            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        AmrIntVect_t iv(i, j, k);
                        if ( cells.find(iv) != cells.end() && cells[iv] <= maxNumPart_m )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::tagForMinNumParticles_m(int lev, TagBoxArray_t& tags,
                                        AmrReal_t time, int ngrow)
{
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // we need to update first
    amrpbase_p->update(0, lev, true);
    const auto& LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();
    
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);
    
    // count particles per cell
    std::map<AmrIntVect_t, size_t> cells;
    for (size_t i = lBegin; i < lEnd; ++i) {
        AmrIntVect_t iv = layout_mp->Index(bunch_mp->R[i], lev);
        ++cells[iv];
    }
        
    const int clearval = TagBox_t::CLEAR;
    const int   tagval = TagBox_t::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter_t mfi(*rho_m[lev], true); mfi.isValid(); ++mfi) {
            
            const Box_t&  tilebx  = mfi.tilebox();
            TagBox_t&     tagfab  = tags[mfi];

            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        AmrIntVect_t iv(i, j, k);
                        if ( cells.find(iv) != cells.end() &&
                            cells[iv] >= minNumPart_m )
                        {
                            tagfab(iv) = tagval;
                        } else
                            tagfab(iv) = clearval;
                    }
                }
            }
        }
    }
}


void AmrBoxLib::initBaseLevel_m(const AmrIntArray_t& nGridPts) {
    // we need to set the AmrCore variable
    finest_level = 0;
    
    AmrIntVect_t low(0, 0, 0); 
    AmrIntVect_t high(nGridPts[0] - 1,
                      nGridPts[1] - 1,
                      nGridPts[2] - 1);
    
    const Box_t bx(low, high);
    AmrGrid_t ba(bx);
    ba.maxSize( this->maxGridSize(0) );
    
    // chop grids to guarantee #grids == #procs on base level
    ba = this->MakeBaseGrids();

    AmrProcMap_t dmap;
    dmap.define(ba, amrex::ParallelDescriptor::NProcs());
    
    this->RemakeLevel(0, 0.0, ba, dmap);
    
    layout_mp->define(this->geom);
    layout_mp->define(this->ref_ratio);    
}


void AmrBoxLib::initParmParse_m(const AmrInfo& info, AmrLayout_t* layout_p) {
    
    /*
     * All parameters that we set with the OPAL input file
     */
    amrex::ParmParse pAmr("amr");
    pAmr.add("max_grid_size_x", info.maxgrid[0]);
    pAmr.add("max_grid_size_y", info.maxgrid[1]);
    pAmr.add("max_grid_size_z", info.maxgrid[2]);
    
    pAmr.add("blocking_factor_x", info.bf[0]);
    pAmr.add("blocking_factor_y", info.bf[1]);
    pAmr.add("blocking_factor_z", info.bf[2]);
    
    const int nratios_vect = info.maxlevel * AMREX_SPACEDIM;
    
    AmrIntArray_t refratio(nratios_vect);
    
    for (int i = 0; i < info.maxlevel; ++i) {
        refratio[i * AMREX_SPACEDIM]     = info.refratio[0];
        refratio[i * AMREX_SPACEDIM + 1] = info.refratio[1];
        refratio[i * AMREX_SPACEDIM + 2] = info.refratio[2];
    }
    
    pAmr.addarr("ref_ratio_vect", refratio);
    
    amrex::ParmParse pGeom("geometry");
    AmrIntArray_t isPeriodic = {
        layout_p->Geom(0).isPeriodic(0),
        layout_p->Geom(0).isPeriodic(1),
        layout_p->Geom(0).isPeriodic(2)
    };
    pGeom.addarr("is_periodic", isPeriodic);
    
    
    /*
     * "All" other parameters, we moslty take the
     * defaults of the code
     * 
     * ATTENTION Be careful with default values since they might
     *           be changed by the AMReX developers!
     * 
     * Parmparse parameters not to be set because
     * alreday given due to constructor
     * ----------------------------------
     * 
     *  AMReX_AmrMesh:
     *      - amr.max_level
     * 
     *  AMReX_Geometry:
     *      - geometry.coord_sys
     *      - geometry.prob_lo
     *      - geometry.prob_hi
     *      - geometry.is_periodic
     *      - geometry.spherical_origin_fix [default: 0]
     *        (not set because we're using Cartesian coordinates)
     * 
     * 
     * ParmParse left out
     * ------------------
     *  AMReX_FabArrayBase:
     *      - fabarray.mfiter_tile_size [default: IntVect(1024000,8,8)]
     *      - fabarray.mfghostiter_tile_size [default: IntVect(1024000, 8, 8)]
     *      - fabarray.comm_tile_size [default: IntVect(1024000, 8, 8)]
     *      - fabarray.maxcomp [default: 25]
     *      - fabarray.do_async_sends [default: true]
     * 
     *  AMReX_MemPool:
     *      - fab.init_snan [default: 0, if not BL_TESTING or DEBUG enabled]
     * 
     *  AMReX_MemProfiler:
     *      - amrex.memory.log [default: "memlog"]
     * 
     *  AMReX_VisMF:
     *      - vismf.v [verbose, default: 0]
     *      - vismf.headerversion [default: VisMF::Header::Version_v1 = 1]
     *      - vismf.groupsets [default: false]
     *      - vismf.setbuf [default: true]
     *      - vismf.usesingleread [default: false]
     *      - vismf.usesinglewrite [default: false]
     *      - vismf.checkfilepositions [default: false]
     *      - vismf.usepersistentifstreams [default: true]
     *      - vismf.usesynchronousreads [default: false]
     *      - vismf.usedynamicsetselection [default: true]
     *      - vismf.iobuffersize [default: VisMF::IO_Buffer_Size = 262144 * 8]
     * 
     *  AMReX_ParallelDescriptor:
     *      Needs change of ParallelDescriptor::SetNProcsSidecars in the function
     *      ParallelDescriptor::StartParallel()
     *      - amrex.ncolors [default: m_nCommColors = 1]
     *      - team.size     [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     *      - team.reduce   [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     *      - mpi.onesided  [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     * 
     *  AMReX:
     *      Debugging purposes
     *      - amrex.v [verbose, default: 0]
     *      - amrex.verbose [default: 0]
     *      - amrex.fpe_trap_invalid [invalid, default: 0]
     *      - amrex.fpe_trap_zero [divbyzero, default: 0]
     *      - amrex.fpe_trap_overflow [overflow, 0]
     * 
     *  AMReX_FArrayBox:
     *      For output only
     *      - fab.format [default: FABio::FAB_NATIVE]
     *      For reading in an old FAB
     *      - fab.ordering
     *      - fab.initval [default: quiet_NaN() or max()]
     *      - fab.do_initval [default: false, if not DEBUG or BL_TESTING enabled]
     *      - fab.init_snan [default: false, if not DEBUG or BL_TESTING enabled]
     * 
     *  AMReX_MCMultiGrid:
     *      We do NOT use this solver
     *      - mg.maxiter [default: 40]
     *      - mg.numiter [default: -1]
     *      - mg.nu_0 [default: 1]
     *      - mg.nu_1 [default: 2]
     *      - mg.nu_2 [default: 2]
     *      - mg.nu_f [default: 8]
     *      - mg.v [verbose, default: 0]
     *      - mg.usecg [default: 1]
     *      - mg.rtol_b [default: 0.01]
     *      - mg.bot_atol [default: -1.0]
     *      - mg.nu_b [default: 0]
     *      - mg.numLevelsMAX [default: 1024]
     * 
     *  AMReX_MCLinOp:
     *      We do NOT use this solver
     *      - MCLp.harmavg [default: 0]
     *      - MCLp.v [verbose, default: 0]
     *      - MCLp.maxorder [default: 2]
     * 
     *  AMReX_MCCGSolver:
     *      We do NOT use this solver
     *      - cg.maxiter [default: 40]
     *      - cg.v [verbose, default: 0]
     *      - cg.isExpert [default: 0]
     *      - MCCGSolver::def_unstable_criterion [default: 10]
     * 
     *  AMReX_LinOp:
     *      We do NOT use this class
     *      - Lp.harmavg [default: 0]
     *      - Lp.v [verbose, default: 0]
     *      - Lp.maxorder [default: 2]
     *      - LinOp::LinOp_grow [default: 1]
     * 
     *  AMReX_MultiGrid (single level solver):
     *      We do NOT use this solver, see BoxLibSolvers/FMGPoissonSolver for
     *      explanation of parameters
     *      - mg.v [verbose, default: 0]
     *      - mg.nu_0 [default: 1]
     *      - mg.nu_1 [default: 2]
     *      - mg.nu_2 [default: 2]
     *      - mg.nu_f [default: 8]
     *      - mg.nu_b [default: 0]
     *      - mg.usecg [default: 1]
     *      - mg.rtol_b [default: 0.0001 or 0.01 (if CG_USE_OLD_CONVERGENCE_CRITERIA enabled)]
     *      - mg.verbose [default: 0]
     *      - mg.maxiter [default: 40]
     *      - mg.bot_atol [default: -1.0]
     *      - mg.maxiter_b [default: 120]
     *      - mg.numLevelsMAX [default: 1024]
     *      - mg.smooth_on_cg_unstable [default: 1]
     *      - mg.use_Anorm_for_convergence [default: 1]
     * 
     *  AMReX_CGSolver (single level solver):
     *      We do NOT use this solver
     *      - cg.v [verbose, default: 0]
     *      - cg.SSS [default: SSS_MAX = 4]
     *      - cg.maxiter [default: 80]
     *      - cg.verbose [default: 0]
     *      - cg.variable_SSS [default: true]
     *      - cg.use_jbb_precond [default: 0]
     *      - cg.use_jacobi_precond [default: 0]
     *      - cg.unstable_criterion [default: 10]
     *      - CGSolver::def_cg_solver [default: BiCGStab]
     *      - cg.cg_solver [default: CGSolver::def_cg_solver]
     * 
     *  AMReX_Amr:
     *      We do NOT use this class
     *      - amr.regrid_on_restart [default: 0]
     *      - amr.use_efficient_regrid [default: 0]
     *      - amr.plotfile_on_restart [default: 0]
     *      - amr.checkpoint_on_restart [default: 0]
     *      - amr.compute_new_dt_on_regrid [default: 0]
     *      - amr.mffile_nstreams [default: 1]
     *      - amr.probinit_natonce [default: 32]
     *      - amr.file_name_digits [default: 5]
     *      - amr.initial_grid_file
     *      - amr.regrid_file
     *      - amr.message_int [default: 10]
     *      - amr.run_log [default: not parsed]
     *      - amr.run_log_terse [default: not parsed]
     *      - amr.grid_log [default: not parsed]
     *      - amr.data_log [default: not parsed]
     *      - amr.probin_file [default: not parsed]
     *      - amr.restart
     *      - amr.restart_from_plotfile
     *      - amr.regrid_int [default: 1 at all levels]
     *      - amr.rebalance_grids [default: 0]
     *      - amr.nosub [default: not parsed]
     *      - amr.subcycling_mode [default: "Auto"]
     *      - amr.subcycling_iterations [default: 0]
     *      - amr.checkpoint_files_output
     *      - amr.plot_files_output
     *      - amr.plot_nfiles
     *      - amr.checkpoint_nfiles
     *      - amr.check_file [default: "chk"]
     *      - amr.check_int [default: -1]
     *      - amr.check_per [default: -1.0] 
     *      - amr.plot_file [default: "plt"]
     *      - amr.plot_int [default: -1]
     *      - amr.plot_per [default: -1.0]
     *      - amr.small_plot_file [default: "smallplt"]
     *      - amr.small_plot_int [default: -1]
     *      - amr.small_plot_per [default: -1.0]
     *      - amr.write_plotfile_with_checkpoint [default: 1]
     *      - amr.stream_max_tries [default: 4]
     *      - amr.abort_on_stream_retry_failure [default: false]
     *      - amr.precreateDirectories
     *      - amr.prereadFAHeaders
     *      - amr.plot_headerversion
     *      - amr.checkpoint_headerversion
     * 
     *  AMReX_SlabStat:
     *      We do NOT use this class
     *      - slabstat.boxes
     * 
     *  AMReX_AmrLevel:
     *      We do NOT use this class
     *      - amr.plot_vars [default: not parsed]
     *      - amr.derive_plot_vars [default: not parsed]
     *      - amr.small_plot_vars [default: not parsed]
     * 
     * AMReX_StationData:
     *      We do NOT use this class
     *      - StationData.vars [default: not parsed]
     *      - StationData.coord [default: not parsed]
     */
    
    // verbosity
    pAmr.add("v", 0);
    
    // # cells required for proper nesting.
    pAmr.add("n_proper", 1);
    
    // Grid efficiency. (default: 0.7)
    pAmr.add("grid_eff", 0.95);
    
    int nlev = info.maxlevel + 1;
    
    // Buffer cells around each tagged cell.
    AmrIntArray_t error_buf(nlev, 4 /*buffer*/);
    error_buf[0] = 0;
    pAmr.addarr("n_error_buf", error_buf);
    
    // chop up grids to have more grids than the number of procs
    pAmr.add("refine_grid_layout", true);
    
    
    /*
     * ParmParse for DistributionMapping
     * 
     * round-robin: FAB i is owned by CPU i%N where N is total number of CPUs
     * knapsack:    FABs are partitioned across CPUs such that the total volume
     *              of the Boxes in the underlying BoxArray are as equal across
     *              CPUs as is possible.
     * SFC:         Is based on a space filling curve (default)
     */
    amrex::ParmParse pDmap("DistributionMapping");
    
    // verbosity
    pDmap.add("v", 0);
    
    // max_efficiency   = 0.9
    pDmap.add("efficiency", 0.9);
    
    // SFC = space filling curve
    pDmap.add("sfc_threshold", 0);
    
    /* if > 0:
     *  nworkers = node_size
     *  nteams   = nprocs / node_size
     *  
     *  if nwokers * nteams != nprocs:
     *      nteams   = nprocs
     *      nworkers = 1
     */
    pDmap.add("node_size", 0);
    
    /* Possibilities:
     *  - ROUNDROBIN
     *  - KNAPSACK
     *  - SFC
     *  - PFC
     *  - RRSFC
     */
    pDmap.add("strategy", "SFC");
}


void AmrBoxLib::fillPhysbc_m(AmrField_t& mf, int lev) {
    /* Copied from one of the miniapps:
     * 
     * amrex/Src/AmrTask/tutorials/MiniApps/HeatEquation/physbc.cpp
     * 
     * AMReX - Revision:
     * 
     * commit 8174212e898677d6413cf5e4db44148a52b4b732
     * Author: Weiqun Zhang <weiqunzhang@lbl.gov>
     * Date:   Mon Jul 2 10:40:21 2018 -0700
     */
    if (AmrGeometry_t::isAllPeriodic())
        return;
    // Set up BC; see Src/Base/AMReX_BC_TYPES.H for supported types
    amrex::Vector<amrex::BCRec> bc(mf.nComp());
    for (int n = 0; n < mf.nComp(); ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            if (AmrGeometry_t::isPeriodic(idim))
            {
                bc[n].setLo(idim, amrex::BCType::int_dir); // interior
                bc[n].setHi(idim, amrex::BCType::int_dir);
            }
            else
            {
                bc[n].setLo(idim, amrex::BCType::foextrap); // first-order extrapolation.
                bc[n].setHi(idim, amrex::BCType::foextrap);
            }
        }
    }
    
    amrex::FillDomainBoundary(mf, this->geom[lev], bc);
}
