#ifndef PART_BUNCH_BASE_HPP
#define PART_BUNCH_BASE_HPP

#include "PartBunchBase.h"

#include "Distribution/Distribution.h"

#include "AbstractObjects/OpalData.h"   // OPAL file
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

extern Inform *gmsg;

template <class T, unsigned Dim>
PartBunchBase<T, Dim>::PartBunchBase(AbstractParticle<T, Dim>* pb)
    : R(*(pb->R_p)),
      ID(*(pb->ID_p)),
      myNode_m(Ippl::myNode()),
      nodes_m(Ippl::getNodes()),
      fixed_grid(false),
      pbin_m(nullptr),
      lossDs_m(nullptr),
      pmsg_m(nullptr),
      f_stream(nullptr),
      lowParticleCount_m(false),
//       reference(ref), //FIXME
      unit_state_(units),
      stateOfLastBoundP_(unitless),
      moments_m(),
      dt_m(0.0),
      t_m(0.0),
      eKin_m(0.0),
      dE_m(0.0),
      spos_m(0.0),
      globalMeanR_m(Vector_t(0.0, 0.0, 0.0)),
      globalToLocalQuaternion_m(Quaternion_t(1.0, 0.0, 0.0, 0.0)),
      rmax_m(0.0),
      rmin_m(0.0),
      rrms_m(0.0),
      prms_m(0.0),
      rmean_m(0.0),
      pmean_m(0.0),
      eps_m(0.0),
      eps_norm_m(0.0),
      halo_m(Vector_t(0.0, 0.0, 0.0)),
      rprms_m(0.0),
      Dx_m(0.0),
      Dy_m(0.0),
      DDx_m(0.0),
      DDy_m(0.0),
      hr_m(-1.0),
      nr_m(0),
      fs_m(nullptr),
      couplingConstant_m(0.0),
      qi_m(0.0),
      distDump_m(0),
      fieldDBGStep_m(0),
      dh_m(1e-12),
      tEmission_m(0.0),
      bingamma_m(nullptr),
      binemitted_m(nullptr),
      lPath_m(0.0),
      stepsPerTurn_m(0),
      localTrackStep_m(0),
      globalTrackStep_m(0),
      numBunch_m(1),
      SteptoLastInj_m(0),
      globalPartPerNode_m(nullptr),
      dist_m(nullptr),
      dcBeam_m(false),
      periodLength_m(Physics::c / 1e9),
      pbase(pb)
{
    setup(pb);

    boundpTimer_m = IpplTimings::getTimer("Boundingbox");
    boundpBoundsTimer_m = IpplTimings::getTimer("Boundingbox-bounds");
    boundpUpdateTimer_m = IpplTimings::getTimer("Boundingbox-update");

    statParamTimer_m = IpplTimings::getTimer("Compute Statistics");
    selfFieldTimer_m = IpplTimings::getTimer("SelfField total");

    histoTimer_m = IpplTimings::getTimer("Histogram");

    distrCreate_m = IpplTimings::getTimer("Create Distr");
    distrReload_m = IpplTimings::getTimer("Load Distr");

    globalPartPerNode_m = std::unique_ptr<size_t[]>(new size_t[Ippl::getNodes()]);

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(std::string("GlobalLosses"), !Options::asciidump));

    pmsg_m.release();
    //    f_stream.release();
    /*
      if(Ippl::getNodes() == 1) {
          f_stream = std::unique_ptr<ofstream>(new ofstream);
          f_stream->open("data/dist.dat", ios::out);
          pmsg_m = std::unique_ptr<Inform>(new Inform(0, *f_stream, 0));
      }
    */
}

template <class T, unsigned Dim>
PartBunchBase<T, Dim>::PartBunchBase(AbstractParticle<T, Dim>* pb, const PartData *ref)
    : PartBunchBase(pb)
{
    reference = ref;
}

/*
 * Bunch common member functions
 */

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::getIfBeamEmitting() {
    if (dist_m != NULL) {
        size_t isBeamEmitted = dist_m->getIfDistEmitting();
        reduce(isBeamEmitted, isBeamEmitted, OpAddAssign());
        if (isBeamEmitted > 0)
            return true;
        else
            return false;
    } else
        return false;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getLastEmittedEnergyBin() {
    /*
     * Get maximum last emitted energy bin.
     */
    int lastEmittedBin = dist_m->getLastEmittedEnergyBin();
    reduce(lastEmittedBin, lastEmittedBin, OpMaxAssign());
    return lastEmittedBin;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getNumberOfEmissionSteps() {
    return dist_m->getNumberOfEmissionSteps();
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getNumberOfEnergyBins() {
    return dist_m->getNumberOfEnergyBins();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::Rebin() {

    size_t isBeamEmitting = dist_m->getIfDistEmitting();
    reduce(isBeamEmitting, isBeamEmitting, OpAddAssign());
    if (isBeamEmitting > 0) {
        *gmsg << "*****************************************************" << endl
              << "Warning: attempted to rebin, but not all distribution" << endl
              << "particles have been emitted. Rebin failed." << endl
              << "*****************************************************" << endl;
    } else {
        if (dist_m->Rebin())
            this->Bin = 0;
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setEnergyBins(int numberOfEnergyBins) {
    bingamma_m = std::unique_ptr<double[]>(new double[numberOfEnergyBins]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[numberOfEnergyBins]);
    for(int i = 0; i < numberOfEnergyBins; i++)
        binemitted_m[i] = 0;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveEnergyBins() {

    if (dist_m != NULL)
        return dist_m->getNumberOfEnergyBins() > 0;
    else
        return false;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::switchToUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == unitless)
        throw SwitcherError("PartBunch::switchToUnitlessPositions",
                            "Cannot make a unitless PartBunch unitless");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] /= Vector_t(Physics::c * dt);
    }

    unit_state_ = unitless;

    if(hasToReset) R.resetDirtyFlag();
}


//FIXME: unify methods, use convention that all particles have own dt
template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::switchOffUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == units)
        throw SwitcherError("PartBunch::switchOffUnitlessPositions",
                            "Cannot apply units twice to PartBunch");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] *= Vector_t(Physics::c * dt);
    }

    unit_state_ = units;

    if(hasToReset) R.resetDirtyFlag();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::do_binaryRepart() {
    // do nothing here
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setDistribution(Distribution *d,
                                            std::vector<Distribution *> addedDistributions,
                                            size_t &np)
{
    Inform m("setDistribution " );
    dist_m = d;

    dist_m->createOpalT(this, addedDistributions, np);

//    if (Options::cZero)
//        dist_m->create(this, addedDistributions, np / 2);
//    else
//        dist_m->create(this, addedDistributions, np);
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::isGridFixed() {
    return fixed_grid;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::hasBinning() {
    return (pbin_m != nullptr);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setTEmission(double t) {
    tEmission_m = t;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getTEmission() {
    return tEmission_m;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::doEmission() {
    if (dist_m != NULL)
        return dist_m->getIfDistEmitting();
    else
        return false;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveBins() const {

    if(pbin_m != NULL)
        return pbin_m->weHaveBins();
    else
        return false;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBins *pbin) {
    pbin_m = pbin;
    *gmsg << *pbin_m << endl;
    setEnergyBins(pbin_m->getNBins());
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBinsCyc *pbin) {
    pbin_m = pbin;
    setEnergyBins(pbin_m->getNBins());
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::emitParticles(double eZ) {
    return dist_m->emitParticles(this, eZ);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::updateNumTotal() {
    size_t numParticles = getLocalNum();
    reduce(numParticles, numParticles, OpAddAssign());
    setTotalNum(numParticles);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::rebin() {
    this->Bin = 0;
    pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = NULL;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getNumBins() {
    if(pbin_m != NULL)
        return pbin_m->getNBins();
    else
        return 0;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getLastemittedBin() {
    if(pbin_m != NULL)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalBinCount(size_t num, int bin) {
    if(pbin_m != NULL) {
        pbin_m->setPartNum(bin, num);
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas() {

    const int emittedBins = dist_m->getNumberOfEnergyBins();

    size_t s = 0;

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;

    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));

    std::unique_ptr<size_t[]> particlesInBin(new size_t[emittedBins]);
    reduce(bingamma_m.get(), bingamma_m.get() + emittedBins, bingamma_m.get(), OpAddAssign());
    reduce(binemitted_m.get(), binemitted_m.get() + emittedBins, particlesInBin.get(), OpAddAssign());
    for(int i = 0; i < emittedBins; i++) {
        size_t &pInBin = particlesInBin[i];
        if(pInBin != 0) {
            bingamma_m[i] /= pInBin;
            INFOMSG(level2 << "Bin " << std::setw(3) << i
                           << " gamma = " << std::setw(8) << std::scientific
                           << std::setprecision(5) << bingamma_m[i]
                           << "; NpInBin= " << std::setw(8)
                           << std::setfill(' ') << pInBin << endl);
        } else {
            bingamma_m[i] = 1.0;
            INFOMSG(level2 << "Bin " << std::setw(3) << i << " has no particles " << endl);
        }
        s += pInBin;
    }
    particlesInBin.reset();


    if(s != getTotalNum() && !OpalData::getInstance()->hasGlobalGeometry())
        ERRORMSG("sum(Bins)= " << s << " != sum(R)= " << getTotalNum() << endl;);

    if(emittedBins >= 2) {
        for(int i = 1; i < emittedBins; i++) {
            if(binemitted_m[i - 1] != 0 && binemitted_m[i] != 0)
                INFOMSG(level2 << "d(gamma)= " << 100.0 * std::abs(bingamma_m[i - 1] - bingamma_m[i]) / bingamma_m[i] << " [%] "
                        << "between bin " << i - 1 << " and " << i << endl);
        }
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas_cycl() {

    const int emittedBins = pbin_m->getLastemittedBin();

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;
    for(unsigned int n = 0; n < getLocalNum(); n++) {
        if ( this->Bin[n] > -1 ) {
            bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));
	}
    }

    allreduce(*bingamma_m.get(), emittedBins, std::plus<double>());

    for(int i = 0; i < emittedBins; i++) {
        if(pbin_m->getTotalNumPerBin(i) > 0)
            bingamma_m[i] /= pbin_m->getTotalNumPerBin(i);
        else
            bingamma_m[i] = 0.0;
        INFOMSG("Bin " << i << " : particle number = " << pbin_m->getTotalNumPerBin(i)
                       << " gamma = " << bingamma_m[i] << endl);
    }

}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getBinGamma(int bin) {
    return bingamma_m[bin];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBinCharge(int bin, double q) {
  this->Q = where(eq(this->Bin, bin), q, 0.0);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBinCharge(int bin) {
  this->Q = where(eq(this->Bin, bin), this->qi_m, 0.0);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::calcNumPartsOutside(Vector_t x) {

    std::size_t localnum = 0;
    const Vector_t meanR = get_rmean();

    for(unsigned long k = 0; k < getLocalNum(); ++ k)
        if (std::abs(R[k](0) - meanR(0)) > x(0) ||
            std::abs(R[k](1) - meanR(1)) > x(1) ||
            std::abs(R[k](2) - meanR(2)) > x(2)) {

            ++localnum;
        }

    gather(&localnum, &globalPartPerNode_m[0], 1);

    size_t npOutside = std::accumulate(globalPartPerNode_m.get(),
                                       globalPartPerNode_m.get() + Ippl::getNodes(), 0,
                                       std::plus<size_t>());

    return npOutside;
}


/**
 * \method calcLineDensity()
 * \brief calculates the 1d line density (not normalized) and append it to a file.
 * \see ParallelTTracker
 * \warning none yet
 *
 * DETAILED TODO
 *
 */
template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcLineDensity(unsigned int nBins,
                                            std::vector<double> &lineDensity,
                                            std::pair<double, double> &meshInfo)
{
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    if (nBins < 2) {
        Vektor<int, 3>/*NDIndex<3>*/ grid;
        this->updateDomainLength(grid);
        nBins = grid[2];
    }

    double length = rmax(2) - rmin(2);
    double zmin = rmin(2) - dh_m * length, zmax = rmax(2) + dh_m * length;
    double hz = (zmax - zmin) / (nBins - 2);
    double perMeter = 1.0 / hz;//(zmax - zmin);
    zmin -= hz;

    lineDensity.resize(nBins, 0.0);
    std::fill(lineDensity.begin(), lineDensity.end(), 0.0);

    const unsigned int lN = getLocalNum();
    for (unsigned int i = 0; i < lN; ++ i) {
        const double z = R[i](2) - 0.5 * hz;
        unsigned int idx = (z - zmin) / hz;
        double tau = (z - zmin) / hz - idx;

        lineDensity[idx] += Q[i] * (1.0 - tau) * perMeter;
        lineDensity[idx + 1] += Q[i] * tau * perMeter;
    }

    reduce(&(lineDensity[0]), &(lineDensity[0]) + nBins, &(lineDensity[0]), OpAddAssign());

    meshInfo.first = zmin;
    meshInfo.second = hz;
}



template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::boundp() {
    /*
      Assume rmin_m < 0.0
     */

    IpplTimings::startTimer(boundpTimer_m);
    //if(!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;

    if(!isGridFixed()) {
        const int dimIdx = (dcBeam_m? 2: 3);

        /**
            In case of dcBeam_m && hr_m < 0
            this is the first call to boundp and we
            have to set hr completely i.e. x,y and z.
         */

        this->updateDomainLength(nr_m);
        IpplTimings::startTimer(boundpBoundsTimer_m);
        get_bounds(rmin_m, rmax_m);
        IpplTimings::stopTimer(boundpBoundsTimer_m);
        Vector_t len = rmax_m - rmin_m;

        double volume = 1.0;
        for(int i = 0; i < dimIdx; i++) {
            double length = std::abs(rmax_m[i] - rmin_m[i]);
            if (length < 1e-10) {
                rmax_m[i] += 1e-10;
                rmin_m[i] -= 1e-10;
            } else {
                rmax_m[i] += dh_m * length;
                rmin_m[i] -= dh_m * length;
            }
            hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
        }
        if (dcBeam_m) {
            rmax_m[2] = rmin_m[2] + periodLength_m;
            hr_m[2] = periodLength_m / (nr_m[2] - 1);
        }
        for (int i = 0; i < dimIdx; ++ i) {
            volume *= std::abs(rmax_m[i] - rmin_m[i]);
        }

        if (getIfBeamEmitting() && dist_m != NULL) {
            // keep particles per cell ratio high, don't spread a hand full particles across the whole grid
            double percent = std::max(1.0 / (nr_m[2] - 1), dist_m->getPercentageEmitted());
            double length  = std::abs(rmax_m[2] - rmin_m[2]) / (1.0 + 2 * dh_m);
            if (percent < 1.0 && percent > 0.0) {
                rmax_m[2] -= dh_m * length;
                rmin_m[2] = rmax_m[2] - length / percent;

                length /= percent;

                rmax_m[2] += dh_m * length;
                rmin_m[2] -= dh_m * length;

                hr_m[2] = (rmax_m[2] - rmin_m[2]) / (nr_m[2] - 1);
            }
        }

        if (volume < 1e-21 && getTotalNum() > 1 && std::abs(sum(Q)) > 0.0) {
            WARNMSG(level1 << "!!! Extremely high particle density detected !!!" << endl);
        }
        //INFOMSG("It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " << rmin_m << endl);

        if(hr_m[0] * hr_m[1] * hr_m[2] <= 0) {
            throw GeneralClassicException("boundp() ", "h<0, can not build a mesh");
        }

        Vector_t origin = rmin_m - Vector_t(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0);
        this->updateFields(hr_m, origin);
    }
    IpplTimings::startTimer(boundpUpdateTimer_m);
    update();
    IpplTimings::stopTimer(boundpUpdateTimer_m);
    R.resetDirtyFlag();

    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::boundp_destroy() {

    Inform gmsgAll("boundp_destroy ", INFORM_ALL_NODES);

    Vector_t len;
    const int dimIdx = 3;
    IpplTimings::startTimer(boundpTimer_m);

    std::unique_ptr<size_t[]> countLost;
    if(weHaveBins()) {
        const int tempN = pbin_m->getLastemittedBin();
        countLost = std::unique_ptr<size_t[]>(new size_t[tempN]);
        for(int ii = 0; ii < tempN; ii++) countLost[ii] = 0;
    }

    this->updateDomainLength(nr_m);


    IpplTimings::startTimer(boundpBoundsTimer_m);
    get_bounds(rmin_m, rmax_m);
    IpplTimings::stopTimer(boundpBoundsTimer_m);

    len = rmax_m - rmin_m;

    calcBeamParameters();

    int checkfactor = Options::remotePartDel;
    if (checkfactor != 0) {
        //INFOMSG("checkfactor= " << checkfactor << endl);
        // check the bunch if its full size is larger than checkfactor times of its rms size
        if(checkfactor < 0) {
            checkfactor *= -1;
            if (len[0] > checkfactor * rrms_m[0] ||
                len[1] > checkfactor * rrms_m[1] ||
                len[2] > checkfactor * rrms_m[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the distance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (std::abs(R[ii](0) - rmean_m(0)) > checkfactor * rrms_m[0] ||
                        std::abs(R[ii](1) - rmean_m(1)) > checkfactor * rrms_m[1] ||
                        std::abs(R[ii](2) - rmean_m(2)) > checkfactor * rrms_m[2])
                    {
                        // put particle onto deletion list
                        destroy(1, ii);
                        //update bin parameter
                        if (weHaveBins())
                            countLost[Bin[ii]] += 1 ;
                        /* INFOMSG("REMOTE PARTICLE DELETION: ID = " << ID[ii] << ", R = " << R[ii]
                         * << ", beam rms = " << rrms_m << endl;);
                         */
                    }
                }
            }
        }
        else {
            if (len[0] > checkfactor * rrms_m[0] ||
                len[2] > checkfactor * rrms_m[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the distance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (std::abs(R[ii](0) - rmean_m(0)) > checkfactor * rrms_m[0] ||
                        std::abs(R[ii](2) - rmean_m(2)) > checkfactor * rrms_m[2])
                    {
                        // put particle onto deletion list
                        destroy(1, ii);
                        //update bin parameter
                        if (weHaveBins())
                            countLost[Bin[ii]] += 1 ;
                        /* INFOMSG("REMOTE PARTICLE DELETION: ID = " << ID[ii] << ", R = " << R[ii]
                         * << ", beam rms = " << rrms_m << endl;);
                         */
                    }
                }
            }
        }
    }

    for(int i = 0; i < dimIdx; i++) {
        double length = std::abs(rmax_m[i] - rmin_m[i]);
        rmax_m[i] += dh_m * length;
        rmin_m[i] -= dh_m * length;
        hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
    }

    // rescale mesh
    this->updateFields(hr_m, rmin_m);

    if(weHaveBins()) {
        pbin_m->updatePartInBin_cyc(countLost.get());
    }

    IpplTimings::startTimer(boundpUpdateTimer_m);
    update();
    IpplTimings::stopTimer(boundpUpdateTimer_m);

    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::boundp_destroyT() {

    const unsigned int minNumParticlesPerCore = getMinimumNumberOfParticlesPerCore();

    this->updateDomainLength(nr_m);

    std::unique_ptr<size_t[]> tmpbinemitted;

    boundp();

    size_t ne = 0;
    const size_t localNum = getLocalNum();

    if(weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for(int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                ne++;
                destroy(1, i);
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        for(unsigned int i = 0; i < localNum; i++) {
            if((Bin[i] < 0) && ((localNum - ne) > minNumParticlesPerCore)) {   // need in minimum x particles per node
                ne++;
                destroy(1, i);
            }
        }
        lowParticleCount_m = ((localNum - ne) <= minNumParticlesPerCore);
        reduce(lowParticleCount_m, lowParticleCount_m, OpOr());
    }

    if (lowParticleCount_m) {
        Inform m ("boundp_destroyT a) ", INFORM_ALL_NODES);
        m << level3 << "Warning low number of particles on some cores localNum= "
          << localNum << " ne= " << ne << " NLocal= " << getLocalNum() << endl;
    } else {
        boundp();
    }
    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if(weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    reduce(ne, ne, OpAddAssign());
    return ne;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::destroyT() {

    const unsigned int minNumParticlesPerCore = getMinimumNumberOfParticlesPerCore();
    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = getLocalNum();
    const size_t totalNum = getTotalNum();
    size_t ne = 0;

    if(weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for(int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                destroy(1, i);
                ++ ne;
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        for(size_t i = 0; i < localNum; i++) {
            if((Bin[i] < 0)) {
                if ((localNum - ne) > minNumParticlesPerCore) {   // need in minimum x particles per node
                    ne++;
                    destroy(1, i);
                }
            }
        }
        lowParticleCount_m = ((localNum - ne) <= minNumParticlesPerCore);
        reduce(lowParticleCount_m, lowParticleCount_m, OpOr());
    }

    if (ne > 0) {
        performDestroy(true);
    }

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if (weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    size_t newTotalNum = getLocalNum();
    reduce(newTotalNum, newTotalNum, OpAddAssign());

    setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPx(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPy(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPz(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPx0(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPy0(int i) {
    return 0;
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getX(int i) {
    return this->R[i](0);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getY(int i) {
    return this->R[i](1);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getZ(int i) {
    return this->R[i](2);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getX0(int i) {
    return 0.0;
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getY0(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setZ(int i, double zcoo)
{
    // nothing done here
};


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::get_bounds(Vector_t &rmin, Vector_t &rmax) {

    this->getLocalBounds(rmin, rmax);

    double min[2*Dim];

    for (unsigned int i = 0; i < Dim; ++i) {
        min[2*i] = rmin[i];
        min[2*i + 1] = -rmax[i];
    }

    allreduce(min, 2*Dim, std::less<double>());

    for (unsigned int i = 0; i < Dim; ++i) {
        rmin[i] = min[2*i];
        rmax[i] = -min[2*i + 1];
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::getLocalBounds(Vector_t &rmin, Vector_t &rmax) {
    const size_t localNum = getLocalNum();
    if (localNum == 0) {
        double maxValue = 1e8;
        rmin = Vector_t(maxValue, maxValue, maxValue);
        rmax = Vector_t(-maxValue, -maxValue, -maxValue);
	return;
    }

    rmin = R[0];
    rmax = R[0];
    for (size_t i = 1; i < localNum; ++ i) {
        for (unsigned short d = 0; d < 3u; ++ d) {
            if (rmin(d) > R[i](d)) rmin(d) = R[i](d);
	    else if (rmax(d) < R[i](d)) rmax(d) = R[i](d);
        }
    }
}


template <class T, unsigned Dim>
std::pair<Vector_t, double> PartBunchBase<T, Dim>::getBoundingSphere() {
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
std::pair<Vector_t, double> PartBunchBase<T, Dim>::getLocalBoundingSphere() {
    Vector_t rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::push_back(OpalParticle p) {
    Inform msg("PartBunch ");

    create(1);
    size_t i = getTotalNum();

    R[i](0) = p[0];
    R[i](1) = p[1];
    R[i](2) = p[2];

    P[i](0) = p[3];
    P[i](1) = p[4];
    P[i](2) = p[5];

    update();
    msg << "Created one particle i= " << i << endl;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_part(FVector<double, 6> z, int ii) {
    R[ii](0) = z[0];
    P[ii](0) = z[1];
    R[ii](1) = z[2];
    P[ii](1) = z[3];
    R[ii](2) = z[4];
    P[ii](2) = z[5];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_part(OpalParticle p, int ii) {
    R[ii](0) = p[0];
    P[ii](0) = p[1];
    R[ii](1) = p[2];
    P[ii](1) = p[3];
    R[ii](2) = p[4];
    P[ii](2) = p[5];
}


template <class T, unsigned Dim>
OpalParticle PartBunchBase<T, Dim>::get_part(int ii) {
    OpalParticle part;
    part[0] = R[ii](0);
    part[1] = P[ii](0);
    part[2] = R[ii](1);
    part[3] = P[ii](1);
    part[4] = R[ii](2);
    part[5] = P[ii](2);
    return part;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                                              double &axmax, double &aymax) {
    axmax = aymax = 0.0;
    OpalParticle part;

    for(unsigned int ii = 0; ii < getLocalNum(); ii++) {

        part = get_part(ii);

        double xnor =
            D(0, 0) * part.x()  + D(0, 1) * part.px() + D(0, 2) * part.y() +
            D(0, 3) * part.py() + D(0, 4) * part.t()  + D(0, 5) * part.pt();
        double pxnor =
            D(1, 0) * part.x()  + D(1, 1) * part.px() + D(1, 2) * part.y() +
            D(1, 3) * part.py() + D(1, 4) * part.t()  + D(1, 5) * part.pt();
        double ynor =
            D(2, 0) * part.x()  + D(2, 1) * part.px() + D(2, 2) * part.y() +
            D(2, 3) * part.py() + D(2, 4) * part.t()  + D(2, 5) * part.pt();
        double pynor =
            D(3, 0) * part.x()  + D(3, 1) * part.px() + D(3, 2) * part.y() +
            D(3, 3) * part.py() + D(3, 4) * part.t()  + D(3, 5) * part.pt();

        axmax = std::max(axmax, (xnor * xnor + pxnor * pxnor));
        aymax = std::max(aymax, (ynor * ynor + pynor * pynor));
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setdT(double dt) {
    dt_m = dt;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getdT() const {
    return dt_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setT(double t) {
    t_m = t;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getT() const {
    return t_m;
}


template <class T, unsigned Dim>

double PartBunchBase<T, Dim>::get_sPos() {
    if(sum(PType != ParticleType::REGULAR)) {
        const size_t n = getLocalNum();
        size_t numPrimBeamParts = 0;
        double z = 0.0;
        if(n != 0) {
            for(size_t i = 0; i < n; i++) {
                if(PType[i] == ParticleType::REGULAR) {
                    z += R[i](2);
                    numPrimBeamParts++;
                }
            }
        }
        reduce(z, z, OpAddAssign());
	reduce(numPrimBeamParts, numPrimBeamParts, OpAddAssign());
	if(numPrimBeamParts != 0)
            z = z / numPrimBeamParts;
        return z;
    } else {
        return spos_m;
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_sPos(double s) {
    spos_m = s;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_gamma() const {
    return eKin_m / (getM()*1e-6) + 1.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_meanKineticEnergy() const {
    return eKin_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_origin() const {
    return rmin_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_maxExtent() const {
    return rmax_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_centroid() const {
    return rmean_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rrms() const {
    return rrms_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rprms() const {
    return rprms_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rmean() const {
    return rmean_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_prms() const {
    return prms_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_pmean() const {
    return pmean_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_pmean_Distribution() const {
    if (dist_m && dist_m->getType() != DistrTypeT::FROMFILE)
        return dist_m->get_pmean();

    double gamma = 0.1 / getM() + 1; // set default 0.1 eV
    return Vector_t(0, 0, sqrt(std::pow(gamma, 2) - 1));
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_emit() const {
    return eps_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_norm_emit() const {
    return eps_norm_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_halo() const {
    return halo_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_hr() const {
    return hr_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_Dx() const {
    return Dx_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_Dy() const {
    return Dy_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_DDx() const {
    return DDx_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_DDy() const {
    return DDy_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_meshEnlargement(double dh) {
    dh_m = dh;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::gatherLoadBalanceStatistics() {

    for(int i = 0; i < Ippl::getNodes(); i++)
        globalPartPerNode_m[i] = 0;

    std::size_t localnum = getLocalNum();
    gather(&localnum, &globalPartPerNode_m[0], 1);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getLoadBalance(int p) const {
    return globalPartPerNode_m[p];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::get_PBounds(Vector_t &min, Vector_t &max) const {
    bounds(this->P, min, max);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcBeamParameters() {
    using Physics::c;

    Vector_t eps2, fac, rsqsum, psqsum, rpsum;

    IpplTimings::startTimer(statParamTimer_m);

    const size_t locNp = getLocalNum();
    const size_t totalNum = getTotalNum();
    const double zero = 0.0;

    get_bounds(rmin_m, rmax_m);

    if(totalNum == 0) {
        for(unsigned int i = 0 ; i < Dim; i++) {
            rmean_m(i) = 0.0;
            pmean_m(i) = 0.0;
            rrms_m(i) = 0.0;
            prms_m(i) = 0.0;
            eps_norm_m(i)  = 0.0;
        }
        rprms_m = 0.0;
        eKin_m = 0.0;
        eps_m = 0.0;
        IpplTimings::stopTimer(statParamTimer_m);
        return;
    }

    const size_t intN = calcMoments();
    const double N = static_cast<double>(intN);

    for(unsigned int i = 0 ; i < Dim; i++) {
        rmean_m(i) = centroid_m[2 * i] / N;
        pmean_m(i) = centroid_m[(2 * i) + 1] / N;
        rsqsum(i) = moments_m(2 * i, 2 * i) - N * rmean_m(i) * rmean_m(i);
        psqsum(i) = moments_m((2 * i) + 1, (2 * i) + 1) - N * pmean_m(i) * pmean_m(i);
        if(psqsum(i) < 0)
            psqsum(i) = 0;
        rpsum(i) = moments_m((2 * i), (2 * i) + 1) - N * rmean_m(i) * pmean_m(i);
    }
    eps2 = (rsqsum * psqsum - rpsum * rpsum) / (N * N);
    rpsum /= N;

    for(unsigned int i = 0 ; i < Dim; i++) {
        rrms_m(i) = sqrt(rsqsum(i) / N);
        prms_m(i) = sqrt(psqsum(i) / N);
        eps_norm_m(i)  =  std::sqrt(std::max(eps2(i), zero));
        double tmp = rrms_m(i) * prms_m(i);
        fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
    }

    rprms_m = rpsum * fac;

    Dx_m = moments_m(0, 5) / N;
    DDx_m = moments_m(1, 5) / N;

    Dy_m = moments_m(2, 5) / N;
    DDy_m = moments_m(3, 5) / N;

    // Find unnormalized emittance.
    double gamma = 0.0;
    for(size_t i = 0; i < locNp; i++)
        gamma += sqrt(1.0 + dot(P[i], P[i]));

    allreduce(&gamma, 1, std::plus<double>());
    gamma /= N;

    calcEMean();

    // The computation of the energy spread is an estimation
    // based on the standard deviation of the longitudinal
    // momentum:
    // Var[f(P)] ~= (df/dP)(E[P])^2 Var[P]
    const double m0 = getM() * 1.E-6;
    double tmp = 1.0 / std::pow(eKin_m / m0 + 1., 2.0);
    if (OpalData::getInstance()->isInOPALCyclMode()) {
        dE_m = prms_m(1) * m0 * sqrt(1.0 - tmp);
    } else {
        dE_m = prms_m(2) * m0 * sqrt(1.0 - tmp);
    }

    eps_m = eps_norm_m / Vector_t(gamma * sqrt(1.0 - 1.0 / (gamma * gamma)));
    IpplTimings::stopTimer(statParamTimer_m);

}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcBeamParametersInitial() {
    using Physics::c;

    const double N =  static_cast<double>(getTotalNum());

    if(N == 0) {
        rmean_m = Vector_t(0.0);
        pmean_m = Vector_t(0.0);
        rrms_m  = Vector_t(0.0);
        prms_m  = Vector_t(0.0);
        eps_m   = Vector_t(0.0);
    } else {
        if(Ippl::myNode() == 0) {
            // fixme:  the following code is crap!
            // Only use one node as this function will get called only once before
            // particles have been emitted (at least in principle).
            Vector_t eps2, fac, rsqsum, psqsum, rpsum;

            const double zero = 0.0;
            const double  N =  static_cast<double>(pbin_m->getNp());
            calcMomentsInitial();

            for(unsigned int i = 0 ; i < Dim; i++) {
                rmean_m(i) = centroid_m[2 * i] / N;
                pmean_m(i) = centroid_m[(2 * i) + 1] / N;
                rsqsum(i) = moments_m(2 * i, 2 * i) - N * rmean_m(i) * rmean_m(i);
                psqsum(i) = moments_m((2 * i) + 1, (2 * i) + 1) - N * pmean_m(i) * pmean_m(i);
                if(psqsum(i) < 0)
                    psqsum(i) = 0;
                rpsum(i) =  moments_m((2 * i), (2 * i) + 1) - N * rmean_m(i) * pmean_m(i);
            }
            eps2 = (rsqsum * psqsum - rpsum * rpsum) / (N * N);
            rpsum /= N;

            for(unsigned int i = 0 ; i < Dim; i++) {

                rrms_m(i) = sqrt(rsqsum(i) / N);
                prms_m(i) = sqrt(psqsum(i) / N);
                eps_m(i)  = sqrt(std::max(eps2(i), zero));
                double tmp = rrms_m(i) * prms_m(i);
                fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
            }
            rprms_m = rpsum * fac;
        }
    }
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::calcBunchBeamParameters(MultiBunchDump::beaminfo_t& binfo, short bunch) {
    if ( !OpalData::getInstance()->isInOPALCyclMode() ) {
        return false;
    }

    const unsigned long localNum = getLocalNum();

    long int bunchTotalNum = 0;
    unsigned long bunchLocalNum = 0;

    /* container:
     *
     * ekin, <x>, <y>, <z>, <p_x>, <p_y>, <p_z>,
     * <x^2>, <y^2>, <z^2>, <p_x^2>, <p_y^2>, <p_z^2>,
     * <xp_x>, <y_py>, <zp_z>,
     * <x^3>, <y^3>, <z^3>, <x^4>, <y^4>, <z^4>
     */
    std::vector<double> local(7 * Dim + 1);

    for(unsigned long k = 0; k < localNum; ++ k) {
        if ( bunchNum[k] != bunch || ID[k] == 0 ) {
            continue;
        }

        ++bunchTotalNum;
        ++bunchLocalNum;

        // ekin
        local[0] += std::sqrt(dot(P[k], P[k]) + 1.0);

        for (unsigned int i = 0; i < Dim; ++i) {

            double r = R[k](i);

            // <x>, <y>, <z>
            local[i + 1] += r;

            // <p_x>, <p_y, <p_z>
            local[i + Dim + 1] += P[k](i);

            // <x^2>, <y^2>, <z^2>
            double r2 = r * r;
            local[i + 2 * Dim + 1] += r2;

            // <p_x^2>, <p_y^2>, <p_z^2>
            local[i + 3 * Dim + 1] += P[k](i) * P[k](i);

            // <xp_x>, <y_py>, <zp_z>
            local[i + 4 * Dim + 1] += r * P[k](i);

            // <x^3>, <y^3>, <z^3>
            local[i + 5 * Dim + 1] += r2 * r;

            // <x^4>, <y^4>, <z^4>
            local[i + 6 * Dim + 1] += r2 * r2;
        }
    }

    // inefficient
    allreduce(bunchTotalNum, 1, std::plus<long int>());

    if ( bunchTotalNum == 0 )
        return false;

    // ekin
    const double m0 = getM() * 1.0e-6;
    local[0] -= bunchLocalNum;
    local[0] *= m0;

    allreduce(local.data(), local.size(), std::plus<double>());

    double invN = 1.0 / double(bunchTotalNum);
    binfo.ekin = local[0] * invN;

    binfo.time       = getT() * 1e9;  // ns
    binfo.nParticles = bunchTotalNum;

    for (unsigned int i = 0; i < Dim; ++i) {

        double w = local[i + 1] * invN;
        double pw = local[i + Dim + 1] * invN;
        double w2 = local[i + 2 * Dim + 1] * invN;
        double pw2 = local[i + 3 * Dim + 1] * invN;
        double wpw = local[i + 4 * Dim + 1] * invN;
        double w3 = local[i + 5 * Dim + 1] * invN;
        double w4 = local[i + 6 * Dim + 1] * invN;

        // <x>, <y>, <z>
        binfo.mean[i] = w;

        // central: <p_w^2> - <p_w>^2 (w = x, y, z)
        binfo.prms[i] = pw2 - pw * pw;
        if ( binfo.prms[i] < 0 ) {
            binfo.prms[i] = 0.0;
        }

        // central: <wp_w> - <w><p_w>
        wpw = wpw - w * pw;

        // central: <w^2> - <w>^2 (w = x, y, z)
        binfo.rrms[i] = w2 - w * w;

        // central: normalized emittance
        binfo.emit[i] = (binfo.rrms[i] * binfo.prms[i] - wpw * wpw);
        binfo.emit[i] =  std::sqrt(std::max(binfo.emit[i], 0.0));

        // central: <w^4> - 4 * <w> * <w^3> + 6 * <w>^2 * <w^2> - 3 * <w>^4
        double tmp = w4
                   - 4.0 * w * w3
                   + 6.0 * w * w * w2
                   - 3.0 * w * w * w * w;
        binfo.halo[i] = tmp / ( binfo.rrms[i] * binfo.rrms[i] );

        // central: sqrt(<w^2> - <w>^2) (w = x, y, z)
        binfo.rrms[i] = std::sqrt(binfo.rrms[i]);

        // central: sqrt(<p_w^2> - <p_w>^2)
        binfo.prms[i] = std::sqrt(binfo.prms[i]);
    }

    double tmp = 1.0 / std::pow(binfo.ekin / m0 + 1.0, 2.0);
    binfo.dEkin = binfo.prms[1] * m0 * std::sqrt(1.0 - tmp);

    return true;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getCouplingConstant() const {
    return couplingConstant_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setCouplingConstant(double c) {
    couplingConstant_m = c;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setCharge(double q) {
    qi_m = q;
    if(getTotalNum() != 0)
        Q = qi_m;
    else
        WARNMSG("Could not set total charge in PartBunch::setCharge based on getTotalNum" << endl);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setChargeZeroPart(double q) {
    qi_m = q;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setMass(double mass) {
    M = mass;
}


/// \brief Need Ek for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getEkin() const {
    if(dist_m)
        return dist_m->getEkin();
    else
        return 0.0;
}

/// \brief Need the work function for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getWorkFunctionRf() const {
    if(dist_m)
        return dist_m->getWorkFunctionRf();
    else
        return 0.0;
}
/// \brief Need the laser energy for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getLaserEnergy() const {
    if(dist_m)
        return dist_m->getLaserEnergy();
    else
        return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getCharge() const {
    return sum(Q);
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getChargePerParticle() const {
    return qi_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setSolver(FieldSolver *fs) {
    fs_m = fs;
    fs_m->initSolver(this);

    /**
       CAN not re-inizialize ParticleLayout
       this is an IPPL issue
     */
    if(!OpalData::getInstance()->hasBunchAllocated()) {
        this->initialize(fs_m->getFieldLayout());
//         this->setMesh(fs_m->getMesh());
//         this->setFieldLayout(fs_m->getFieldLayout());
    }
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::hasFieldSolver() {
    if(fs_m)
        return fs_m->hasValidSolver();
    else
        return false;
}


/// \brief Return the fieldsolver type if we have a fieldsolver
template <class T, unsigned Dim>
std::string PartBunchBase<T, Dim>::getFieldSolverType() const {
    if(fs_m)
        return fs_m->getFieldSolverType();
    else
        return "";
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLPath(double s) {
    lPath_m = s;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getLPath() const {
    return lPath_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setStepsPerTurn(int n) {
    stepsPerTurn_m = n;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getStepsPerTurn() const {
    return stepsPerTurn_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setGlobalTrackStep(long long n) {
    globalTrackStep_m = n;
}


template <class T, unsigned Dim>
long long PartBunchBase<T, Dim>::getGlobalTrackStep() const {
    return globalTrackStep_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalTrackStep(long long n) {
    localTrackStep_m = n;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::incTrackSteps() {
    globalTrackStep_m++; localTrackStep_m++;
}


template <class T, unsigned Dim>
long long PartBunchBase<T, Dim>::getLocalTrackStep() const {
    return localTrackStep_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setNumBunch(int n) {
    numBunch_m = n;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getNumBunch() const {
    return numBunch_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setGlobalMeanR(Vector_t globalMeanR) {
    globalMeanR_m = globalMeanR;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::getGlobalMeanR() {
    return globalMeanR_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {

    globalToLocalQuaternion_m = globalToLocalQuaternion;
}


template <class T, unsigned Dim>
Quaternion_t PartBunchBase<T, Dim>::getGlobalToLocalQuaternion() {
    return globalToLocalQuaternion_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setSteptoLastInj(int n) {
    SteptoLastInj_m = n;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getSteptoLastInj() {
    return SteptoLastInj_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calcMeanPhi() {

    const int emittedBins = pbin_m->getLastemittedBin();
    double phi[emittedBins];
    double px[emittedBins];
    double py[emittedBins];
    double meanPhi = 0.0;

    for(int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = 0.0;
        px[ii] = 0.0;
        py[ii] = 0.0;
    }

    for(unsigned int ii = 0; ii < getLocalNum(); ii++) {

        px[Bin[ii]] += P[ii](0);
        py[Bin[ii]] += P[ii](1);
    }

    reduce(px, px + emittedBins, px, OpAddAssign());
    reduce(py, py + emittedBins, py, OpAddAssign());
    for(int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = calculateAngle(px[ii], py[ii]);
        meanPhi += phi[ii];
        INFOMSG("Bin " << ii  << "  mean phi = " << phi[ii] * 180.0 / Physics::pi - 90.0 << "[degree]" << endl);
    }

    meanPhi /= emittedBins;

    INFOMSG("mean phi of all particles " <<  meanPhi * 180.0 / Physics::pi - 90.0 << "[degree]" << endl);

    return meanPhi;
}


// this function reset the BinID for each particles according to its current beta*gamma
// it is for multi-turn extraction cyclotron with small energy gain
// the bin number can be different with the bunch number

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::resetPartBinID2(const double eta) {


    INFOMSG("Before reset Bin: " << endl);
    calcGammas_cycl();
    int maxbin = pbin_m->getNBins();
    size_t partInBin[maxbin];
    for(int ii = 0; ii < maxbin; ii++) partInBin[ii] = 0;

    double pMin0 = 1.0e9;
    double pMin = 0.0;
    double maxbinIndex = 0;

    for(unsigned long int n = 0; n < getLocalNum(); n++) {
        double temp_betagamma = sqrt(pow(P[n](0), 2) + pow(P[n](1), 2));
        if(pMin0 > temp_betagamma)
            pMin0 = temp_betagamma;
    }
    reduce(pMin0, pMin, OpMinAssign());
    INFOMSG("minimal beta*gamma = " << pMin << endl);

    double asinh0 = asinh(pMin);
    for(unsigned long int n = 0; n < getLocalNum(); n++) {

        double temp_betagamma = sqrt(pow(P[n](0), 2) + pow(P[n](1), 2));
#ifdef ENABLE_AMR
        // FIXME: Restrict to 1 GeV proton
        if ( temp_betagamma > 1.808 ) {
            Bin[n] = -1;
            continue;
        }
#endif
        int itsBinID = floor((asinh(temp_betagamma) - asinh0) / eta + 1.0E-6);
        Bin[n] = itsBinID;
        if(maxbinIndex < itsBinID) {
            maxbinIndex = itsBinID;
        }

        if(itsBinID >= maxbin) {
            ERRORMSG("The bin number limit is " << maxbin << ", please increase the energy interval and try again" << endl);
            return false;
        } else
            partInBin[itsBinID]++;

    }

    // partInBin only count particle on the local node.
    pbin_m->resetPartInBin_cyc(partInBin, maxbinIndex);

    // after reset Particle Bin ID, update mass gamma of each bin again
    INFOMSG("After reset Bin: " << endl);
    calcGammas_cycl();

    return true;

}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::resetPartBinBunch() {
    int maxbin = pbin_m->getNBins();
    std::size_t partInBin[maxbin];
    for (int i = 0; i < maxbin; ++i) {
        partInBin[i] = 0;
    }

    for (std::size_t i = 0; i < getLocalNum(); ++i) {
        partInBin[Bin[i]]++;
    }

    // partInBin only count particle on the local node.
    pbin_m->resetPartInBin_cyc(partInBin, numBunch_m - 1);

    // after reset Particle Bin ID, update mass gamma of each bin again
    calcGammas_cycl();

    return true;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getQ() const {
    return reference->getQ();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getM() const {
    return reference->getM();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getP() const {
    return reference->getP();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getE() const {
    return reference->getE();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetQ(double q)  {
    const_cast<PartData *>(reference)->setQ(q);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetM(double m)  {
    const_cast<PartData *>(reference)->setM(m);
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getdE() {
    return dE_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getInitialBeta() const {
    return reference->getBeta();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getInitialGamma() const {
    return reference->getGamma();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getGamma(int i) {
    return 0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getBeta(int i) {
    return 0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::actT()
{
    // do nothing here
};


template <class T, unsigned Dim>
const PartData *PartBunchBase<T, Dim>::getReference() const {
    return reference;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getEmissionDeltaT() {
    return dist_m->getEmissionDeltaT();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::iterateEmittedBin(int binNumber) {
    binemitted_m[binNumber]++;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcEMean() {

    const double totalNp = static_cast<double>(getTotalNum());
    const double locNp = static_cast<double>(getLocalNum());

    eKin_m = 0.0;

    for(unsigned int k = 0; k < locNp; k++) {
        eKin_m += sqrt(dot(P[k], P[k]) + 1.0);
    }

    eKin_m -= locNp;
    eKin_m *= getM() * 1.0e-6;

    reduce(eKin_m, eKin_m, OpAddAssign());

    eKin_m /= totalNp;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::correctEnergy(double avrgp_m) {

    const double totalNp = static_cast<double>(getTotalNum());
    const double locNp = static_cast<double>(getLocalNum());

    double avrgp = 0.0;
    for(unsigned int k = 0; k < locNp; k++)
        avrgp += sqrt(dot(P[k], P[k]));

    reduce(avrgp, avrgp, OpAddAssign());
    avrgp /= totalNp;

    for(unsigned int k = 0; k < locNp; k++)
        P[k](2) =  P[k](2) - avrgp + avrgp_m;
}


template <class T, unsigned Dim>
Inform &PartBunchBase<T, Dim>::print(Inform &os) {

    if(getTotalNum() != 0) {  // to suppress Nan's
        Inform::FmtFlags_t ff = os.flags();

        double lengthUnitConverter = 1;
        double pathLength = get_sPos();
        if (OpalData::getInstance()->isInOPALCyclMode()) {
            lengthUnitConverter = 0.001;
            pathLength = getLPath();
        }

        rmax_m *= lengthUnitConverter;
        rmin_m *= lengthUnitConverter;

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H ********************************************************* \n";
        os << "* NP              = " << getTotalNum() << "\n";
        os << "* Qtot            = " << std::setw(17) << Util::getChargeString(std::abs(sum(Q))) << "         "
        << "Qi    = "             << std::setw(17) << Util::getChargeString(std::abs(qi_m)) << "\n";
        os << "* Ekin            = " << std::setw(17) << Util::getEnergyString(eKin_m) << "         "
           << "dEkin = "             << std::setw(17) << Util::getEnergyString(dE_m) << "\n";
        os << "* rmax            = " << Util::getLengthString(rmax_m, 5) << "\n";
        os << "* rmin            = " << Util::getLengthString(rmin_m, 5) << "\n";
        os << "* rms beam size   = " << Util::getLengthString(rrms_m, 5) << "\n";
        os << "* rms momenta     = " << std::setw(12) << std::setprecision(5) << prms_m << " [beta gamma]\n";
        os << "* mean position   = " << Util::getLengthString(rmean_m, 5) << "\n";
        os << "* mean momenta    = " << std::setw(12) << std::setprecision(5) << pmean_m << " [beta gamma]\n";
        os << "* rms emittance   = " << std::setw(12) << std::setprecision(5) << eps_m << " (not normalized)\n";
        os << "* rms correlation = " << std::setw(12) << std::setprecision(5) << rprms_m << "\n";
        os << "* hr              = " << Util::getLengthString(get_hr(), 5) << "\n";
        os << "* dh              = " << std::setw(13) << std::setprecision(5) << dh_m * 100 << " [%]\n";
        os << "* t               = " << std::setw(17) << Util::getTimeString(getT()) << "         "
           << "dT    = "             << std::setw(17) << Util::getTimeString(getdT()) << "\n";
        os << "* spos            = " << std::setw(17) << Util::getLengthString(pathLength) << "\n";
        os << "* ********************************************************************************** " << endl;
        os.flags(ff);

        rmax_m /= lengthUnitConverter;
        rmin_m /= lengthUnitConverter;
    }
    return os;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::calcMoments() {

    double part[2 * Dim];

    const unsigned long localNum = getLocalNum();

    /* 2 * Dim centroids + Dim * ( 2 * Dim + 1 ) 2nd moments + 2 * Dim (3rd and 4th order moments)
     * --> 1st order moments: 0, ..., 2 * Dim - 1
     * --> 2nd order moments: 2 * Dim, ..., Dim * ( 2 * Dim + 1 )
     * --> 3rd order moments: Dim * ( 2 * Dim + 1 ) + 1, ..., Dim * ( 2 * Dim + 1 ) + Dim
     * (only, <x^3>, <y^3> and <z^3>)
     * --> 4th order moments: Dim * ( 2 * Dim + 1 ) + Dim + 1, ..., Dim * ( 2 * Dim + 1 ) + 2 * Dim
     *
     * For a 6x6 matrix we have each 2nd order moment (except diagonal
     * entries) twice. We only store the upper half of the matrix.
     */
    std::vector<double> loc_moments(4 * Dim + Dim * ( 2 * Dim + 1 ));

    long int totalNum = this->getTotalNum();
    if (OpalData::getInstance()->isInOPALCyclMode()) {

        for(unsigned long k = 0; k < localNum; ++ k) {
            if (ID[k] == 0) {
                part[1] = P[k](0);
                part[3] = P[k](1);
                part[5] = P[k](2);
                part[0] = R[k](0);
                part[2] = R[k](1);
                part[4] = R[k](2);

                unsigned int l = 2 * Dim;
                for (unsigned int i = 0; i < 2 * Dim; ++i) {
                    loc_moments[i] -= part[i];
                    for(unsigned int j = 0; j <= i; j++) {
                        loc_moments[l++] -= part[i] * part[j];
                    }
                }

                for (unsigned int i = 0; i < Dim; ++i) {
                    double r2 = R[k](i) * R[k](i);
                    loc_moments[l] -= r2 * R[k](i);
                    loc_moments[Dim + l++] -= r2 * r2;
                }

                --totalNum;
                break;
            }
        }
        allreduce(totalNum, 1, std::less<long int>());
    }

    for(unsigned long k = 0; k < localNum; ++ k) {
        part[1] = P[k](0);
        part[3] = P[k](1);
        part[5] = P[k](2);
        part[0] = R[k](0);
        part[2] = R[k](1);
        part[4] = R[k](2);

        unsigned int l = 2 * Dim;
        for (unsigned int i = 0; i < 2 * Dim; ++i) {
            loc_moments[i] += part[i];
            for(unsigned int j = 0; j <= i; j++) {
                loc_moments[l++] += part[i] * part[j];
            }
        }

        for (unsigned int i = 0; i < Dim; ++i) {
            double r2 = R[k](i) * R[k](i);
            loc_moments[l] += r2 * R[k](i);
            loc_moments[Dim + l++] += r2 * r2;
        }
    }

    allreduce(loc_moments.data(), loc_moments.size(), std::plus<double>());

    // copy to member variables
    for (unsigned int i = 0; i< 2 * Dim; ++i)
        centroid_m[i] = loc_moments[i];

    unsigned int l = 2 * Dim;
    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        for(unsigned int j = 0; j <= i; j++) {
            moments_m(i, j) = loc_moments[l++];
            moments_m(j, i) = moments_m(i, j);
        }
    }

    /* 4th order central moment: <w^4> - 4<w><w^3> + 6<w>^2<w^2> - 3<w>^4
     * 2nd order central moment: <w^2> - <w>^2
     *
     * with w = x, y, z.
     */
    int j = 2 * Dim + Dim * ( 2 * Dim + 1 );
    double invN = 1.0 / double(totalNum);
    for (unsigned int i = 0; i < Dim; ++i) {
        double w1 = centroid_m[2 * i] * invN;
        double w2 = moments_m(2 * i, 2 * i) * invN;
        double w3 = loc_moments[j + i] * invN;
        double w4 = loc_moments[j + Dim + i] * invN;

        halo_m(i)  =  w4 + w1 * (-4.0 * w3 + w1 * (6.0 * w2 - 3.0 * w1 * w1));
        double tmp = w2 - w1 * w1;
        halo_m(i) /= ( tmp * tmp );
        halo_m(i) -= Options::haloShift;
    }

    return totalNum;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcMomentsInitial() {

    double part[2 * Dim];

    for(unsigned int i = 0; i < 2 * Dim; ++i) {
        centroid_m[i] = 0.0;
        for(unsigned int j = 0; j <= i; ++j) {
            moments_m(i, j) = 0.0;
            moments_m(j, i) = moments_m(i, j);
        }
    }

    for(size_t k = 0; k < pbin_m->getNp(); k++) {
        for(int binNumber = 0; binNumber < pbin_m->getNBins(); binNumber++) {
            std::vector<double> p;

            if(pbin_m->getPart(k, binNumber, p)) {
                part[0] = p.at(0);
                part[1] = p.at(3);
                part[2] = p.at(1);
                part[3] = p.at(4);
                part[4] = p.at(2);
                part[5] = p.at(5);

                for(unsigned int i = 0; i < 2 * Dim; ++i) {
                    centroid_m[i] += part[i];
                    for(unsigned int j = 0; j <= i; ++j) {
                        moments_m(i, j) += part[i] * part[j];
                    }
                }
            }
        }
    }

    for(unsigned int i = 0; i < 2 * Dim; ++i) {
        for(unsigned int j = 0; j < i; ++j) {
            moments_m(j, i) = moments_m(i, j);
        }
    }
}


// angle range [0~2PI) degree
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    return thetaXY >= 0 ? thetaXY : thetaXY + Physics::two_pi;
}


template <class T, unsigned Dim>
Inform &operator<<(Inform &os, PartBunchBase<T, Dim> &p) {
    return p.print(os);
}


/*
 * Virtual member functions
 */

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::runTests() {
    throw OpalException("PartBunchBase<T, Dim>::runTests() ", "No test supported.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetInterpolationCache(bool clearCache) {

}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::swap(unsigned int i, unsigned int j) {
    if (i >= getLocalNum() || j >= getLocalNum() || i == j) return;

    std::swap(R[i], R[j]);
    std::swap(P[i], P[j]);
    std::swap(Q[i], Q[j]);
    std::swap(M[i], M[j]);
    std::swap(Phi[i], Phi[j]);
    std::swap(Ef[i], Ef[j]);
    std::swap(Eftmp[i], Eftmp[j]);
    std::swap(Bf[i], Bf[j]);
    std::swap(Bin[i], Bin[j]);
    std::swap(dt[i], dt[j]);
    std::swap(PType[i], PType[j]);
    std::swap(TriID[i], TriID[j]);
    std::swap(cavityGapCrossed[i], cavityGapCrossed[j]);
    std::swap(bunchNum[i], bunchNum[j]);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCAllPeriodic() {
    throw OpalException("PartBunchBase<T, Dim>::setBCAllPeriodic() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCAllOpen() {
    throw OpalException("PartBunchBase<T, Dim>::setBCAllOpen() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCForDCBeam() {
    throw OpalException("PartBunchBase<T, Dim>::setBCForDCBeam() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::updateFields(const Vector_t& hr, const Vector_t& origin)
{

}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setup(AbstractParticle<T, Dim>* pb) {
    pb->addAttribute(P);
    pb->addAttribute(Q);
    pb->addAttribute(M);
    pb->addAttribute(Phi);
    pb->addAttribute(Ef);
    pb->addAttribute(Eftmp);

    pb->addAttribute(Bf);
    pb->addAttribute(Bin);
    pb->addAttribute(dt);
    pb->addAttribute(PType);
    pb->addAttribute(TriID);
    pb->addAttribute(cavityGapCrossed);

    pb->addAttribute(bunchNum);

    boundpTimer_m       = IpplTimings::getTimer("Boundingbox");
    boundpBoundsTimer_m = IpplTimings::getTimer("Boundingbox-bounds");
    boundpUpdateTimer_m = IpplTimings::getTimer("Boundingbox-update");
    statParamTimer_m    = IpplTimings::getTimer("Compute Statistics");
    selfFieldTimer_m    = IpplTimings::getTimer("SelfField total");

    histoTimer_m        = IpplTimings::getTimer("Histogram");

    distrCreate_m       = IpplTimings::getTimer("Create Distr");
    distrReload_m       = IpplTimings::getTimer("Load Distr");


    globalPartPerNode_m = std::unique_ptr<size_t[]>(new size_t[Ippl::getNodes()]);

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(std::string("GlobalLosses"), !Options::asciidump));

    pmsg_m.release();
    //    f_stream.release();
    /*
      if(Ippl::getNodes() == 1) {
          f_stream = std::unique_ptr<ofstream>(new ofstream);
          f_stream->open("data/dist.dat", ios::out);
          pmsg_m = std::unique_ptr<Inform>(new Inform(0, *f_stream, 0));
      }
    */

    // set the default IPPL behaviour
    setMinimumNumberOfParticlesPerCore(0);
}

template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getTotalNum() const {
    return pbase->getTotalNum();
}

template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getLocalNum() const {
    return pbase->getLocalNum();
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getDestroyNum() const {
    return pbase->getDestroyNum();
}

template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getGhostNum() const {
    return pbase->getGhostNum();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setTotalNum(size_t n) {
    pbase->setTotalNum(n);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalNum(size_t n) {
    pbase->setLocalNum(n);
}

template <class T, unsigned Dim>
unsigned int PartBunchBase<T, Dim>::getMinimumNumberOfParticlesPerCore() const {
    return pbase->getMinimumNumberOfParticlesPerCore();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setMinimumNumberOfParticlesPerCore(unsigned int n) {
    pbase->setMinimumNumberOfParticlesPerCore(n);
}

template <class T, unsigned Dim>
ParticleLayout<T, Dim> & PartBunchBase<T, Dim>::getLayout() {
    return pbase->getLayout();
}

template <class T, unsigned Dim>
const ParticleLayout<T, Dim>& PartBunchBase<T, Dim>::getLayout() const {
    return pbase->getLayout();
}

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::getUpdateFlag(UpdateFlags_t f) const {
    return pbase->getUpdateFlag(f);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setUpdateFlag(UpdateFlags_t f, bool val) {
    pbase->setUpdateFlag(f, val);
}

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::singleInitNode() const {
    return pbase->singleInitNode();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetID() {
    pbase->resetID();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::update() {
    pbase->update();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::update(const ParticleAttrib<char>& canSwap) {
    pbase->update(canSwap);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::createWithID(unsigned id) {
    pbase->createWithID(id);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::create(size_t M) {
    pbase->create(M);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::globalCreate(size_t np) {
    pbase->globalCreate(np);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::destroy(size_t M, size_t I, bool doNow) {
    pbase->destroy(M, I, doNow);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::performDestroy(bool updateLocalNum) {
    pbase->performDestroy(updateLocalNum);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::ghostDestroy(size_t M, size_t I) {
    pbase->ghostDestroy(M, I);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBeamFrequency(double f) {
    periodLength_m = Physics::c / f;
}

template <class T, unsigned Dim>
FMatrix<double, 2 * Dim, 2 * Dim> PartBunchBase<T, Dim>::getSigmaMatrix() {
    //const double  N =  static_cast<double>(this->getTotalNum());

    Vektor<double, 2*Dim> rpmean;
    for (unsigned int i = 0; i < Dim; i++) {
        rpmean(2*i)= rmean_m(i);
        rpmean((2*i)+1)= pmean_m(i);
    }

    FMatrix<double, 2 * Dim, 2 * Dim> sigmaMatrix;// = moments_m / N;
    for (unsigned int i = 0; i < 2 * Dim; i++) {
        for (unsigned int j = 0; j <= i; j++) {
            sigmaMatrix[i][j] = moments_m(i, j) -  rpmean(i) * rpmean(j);
            sigmaMatrix[j][i] = sigmaMatrix[i][j];
        }
    }
    return sigmaMatrix;
}

#endif
