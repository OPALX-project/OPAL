#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

AmrPartBunch::AmrPartBunch(const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), ref)
    , amrobj_mp(nullptr)
    , amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase.get()))
    , fieldlayout_m(nullptr)
{
    amrpbase_mp->initializeAmr();
}

AmrPartBunch::AmrPartBunch(const PartData *ref, pbase_t* pbase_p)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t(&pbase_p->getAmrLayout())), ref)
    , amrobj_mp(nullptr)
    , amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase.get()))
    , fieldlayout_m(nullptr)
{
    amrpbase_mp->initializeAmr();
}

AmrPartBunch::~AmrPartBunch() {
    
}


AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() {
    return amrpbase_mp;
}


const AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() const {
    return amrpbase_mp;
}


void AmrPartBunch::initialize(FieldLayout_t *fLayout) {
//     Layout_t* layout = static_cast<Layout_t*>(&getLayout());
}


void AmrPartBunch::do_binaryRepart() {
    
    if ( amrobj_mp && !amrobj_mp->isRefined() ) {
        /* In the first call to this function we
         * intialize all fine levels
         */
        amrobj_mp->initFineLevels();

    } else {
        bool isForbidTransform = amrpbase_mp->isForbidTransform();

        if ( !isForbidTransform ) {
            /*
             * regrid in boosted frame
             */
            this->updateLorentzFactor();
            // Lorentz transform + mapping to [-1,1]
            amrpbase_mp->domainMapping();
            amrpbase_mp->setForbidTransform(true);
        }

        this->update();

        if ( !isForbidTransform ) {
            amrpbase_mp->setForbidTransform(false);
            // map particles back + undo Lorentz transform
            amrpbase_mp->domainMapping(true);
        }
    }
}


Vector_t AmrPartBunch::get_hr() const {
    const double& scalefactor = amrpbase_mp->getScalingFactor();
    return hr_m * scalefactor;
}


void AmrPartBunch::set_meshEnlargement(double dh) {
    // set dh_m = dh
    PartBunchBase<double, 3>::set_meshEnlargement(dh);
    
    // update base geometry with new mesh enlargement
    AmrLayout_t* layout_p = &amrpbase_mp->getAmrLayout();
    layout_p->setBoundingBox(dh);
    
    // if amrobj_mp != nullptr --> we need to regrid
    this->do_binaryRepart();
}


AmrPartBunch::VectorPair_t AmrPartBunch::getEExtrema() {
    return amrobj_mp->getEExtrema();
}

double AmrPartBunch::getRho(int x, int y, int z) {
    /* This function is called in
     * H5PartWrapperForPC::writeStepData(PartBunchBase<double, 3>* bunch)
     * and
     * H5PartWrapperForPT::writeStepData(PartBunchBase<double, 3>* bunch)
     * in case of Options::rhoDump = true.
     * 
     * Currently, we do not support writing multilevel grid data that's why
     * we average the values down to the coarsest level.
     */
    return amrobj_mp->getRho(x, y, z);
}


FieldLayout_t &AmrPartBunch::getFieldLayout() {
    //TODO Implement
    throw OpalException("AmrPartBunch::getFieldLayout() ", "Not yet Implemented.");
    return *fieldlayout_m;
}


void AmrPartBunch::boundp() {
    IpplTimings::startTimer(boundpTimer_m);
    //if(!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;
    
    if ( amrobj_mp ) {
        /* we do an explicit domain mapping of the particles and then
         * forbid it during the regrid process, this way it's only
         * executed ones --> saves computation
         */
        bool isForbidTransform = amrpbase_mp->isForbidTransform();
            
        if ( !isForbidTransform ) {
            this->updateLorentzFactor();
            // Lorentz transform + mapping to [-1,1]
            amrpbase_mp->domainMapping();
            amrpbase_mp->setForbidTransform(true);
        }
        
        this->update();
        
        if ( !isForbidTransform ) {
            amrpbase_mp->setForbidTransform(false);
            // map particles back + undo Lorentz transform
            amrpbase_mp->domainMapping(true);
        }
        
    } else {
        // At this point an amrobj_mp needs already be set
        throw GeneralClassicException("AmrPartBunch::boundp() ",
                                      "AmrObject pointer is not set.");
    }
    
    R.resetDirtyFlag();
    
    IpplTimings::stopTimer(boundpTimer_m);
}


void AmrPartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    
    if ( !fs_m->hasValidSolver() )
        throw OpalException("AmrPartBunch::computeSelfFields() ",
                            "No field solver.");
    
    amrobj_mp->computeSelfFields();
    
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);
    amrobj_mp->computeSelfFields(bin);
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields_cycl(double gamma) {
    IpplTimings::startTimer(selfFieldTimer_m);
    amrobj_mp->computeSelfFields_cycl(gamma);
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields_cycl(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);
    
    /* make sure it is refined in multi-bunch case
     */
    if ( !amrobj_mp->isRefined() ) {
        amrobj_mp->initFineLevels();
    }
    
    amrobj_mp->computeSelfFields_cycl(bin);
    
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::setAmrDomainRatio(const std::vector<double>& ratio) {
    AmrLayout_t* layout_p = &amrpbase_mp->getAmrLayout();
    layout_p->setDomainRatio(ratio);
}

void AmrPartBunch::gatherLevelStatistics() {
    int nLevel = amrobj_mp->maxLevel() + 1;
    
    std::unique_ptr<size_t[]> partPerLevel( new size_t[nLevel] );
    globalPartPerLevel_m.reset( new size_t[nLevel] );
    
    for (int i = 0; i < nLevel; ++i)
        partPerLevel[i] = globalPartPerLevel_m[i] = 0.0;
    
    // do not modify LocalNumPerLevel in here!!!
    auto& LocalNumPerLevel = amrpbase_mp->getLocalNumPerLevel();
        
    for (size_t i = 0; i < LocalNumPerLevel.size(); ++i)
        partPerLevel[i] = LocalNumPerLevel[i];
    
    reduce(*partPerLevel.get(),
           *globalPartPerLevel_m.get(),
           nLevel, std::plus<size_t>());
}


const size_t& AmrPartBunch::getLevelStatistics(int l) const {
    return globalPartPerLevel_m[l];
}


void AmrPartBunch::updateLorentzFactor(int bin) {
    double gamma = this->get_gamma();

    if ( this->weHaveBins() ) {
        gamma = this->getBinGamma(bin);
    }
    
    
    /* At the beginning of simulation the Lorentz factor
     * is not yet set since PartBunchBase::calcBeamParameters
     * is not yet called. Therefore, we do this ugly workaround.
     */
    bool init = true;
    if ( gamma >= 1.0 ) {
        init = false;
    }
    
    if ( init ) {
        gamma = 1.0;
    }

    updateLorentzFactor(gamma);
}


void AmrPartBunch::updateLorentzFactor(double gamma) {
    // keep all 1.0, except longitudinal direction
    Vector_t lorentzFactor(1.0, 1.0, 1.0);

    if (OpalData::getInstance()->isInOPALCyclMode()) {
        lorentzFactor[1] = gamma;
    } else {
        lorentzFactor[2] = gamma;
    }

    amrpbase_mp->setLorentzFactor(lorentzFactor);
}


void AmrPartBunch::updateFieldContainers_m() {
    
}

void AmrPartBunch::updateDomainLength(Vektor<int, 3>& grid) {
    grid = amrobj_mp->getBaseLevelGridPoints();
}


void AmrPartBunch::updateFields(const Vector_t& hr, const Vector_t& origin) {
    //TODO regrid; called in boundp()
//     amrobj_mp->updateMesh();
}
