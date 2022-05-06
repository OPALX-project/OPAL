//
// Class PartBunchBase
//   Base class for representing particle bunches.
//
// Copyright (c) 2008 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#ifndef PART_BUNCH_BASE_HPP
#define PART_BUNCH_BASE_HPP

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBinsCyc.h"
#include "Algorithms/PartData.h"
#include "Distribution/Distribution.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/SwitcherError.h"
#include "Utilities/Util.h"

#include <cmath>
#include <numeric>

extern Inform* gmsg;

template <class T, unsigned Dim>
PartBunchBase<T, Dim>::PartBunchBase(AbstractParticle<T, Dim>* pb, const PartData* ref)
    : R(*(pb->R_p)),
      ID(*(pb->ID_p)),
      pbin_m(nullptr),
      pmsg_m(nullptr),
      f_stream(nullptr),
      fixed_grid(false),
      reference(ref),
      unit_state_(units),
      stateOfLastBoundP_(unitless),
      moments_m(),
      dt_m(0.0),
      t_m(0.0),
      spos_m(0.0),
      globalMeanR_m(Vector_t(0.0, 0.0, 0.0)),
      globalToLocalQuaternion_m(Quaternion_t(1.0, 0.0, 0.0, 0.0)),
      rmax_m(0.0),
      rmin_m(0.0),
      hr_m(-1.0),
      nr_m(0),
      fs_m(nullptr),
      couplingConstant_m(0.0),
      qi_m(0.0),
      massPerParticle_m(0.0),
      distDump_m(0),
      dh_m(1e-12),
      tEmission_m(0.0),
      bingamma_m(nullptr),
      binemitted_m(nullptr),
      stepsPerTurn_m(0),
      localTrackStep_m(0),
      globalTrackStep_m(0),
      numBunch_m(1),
      bunchTotalNum_m(1),
      bunchLocalNum_m(1),
      SteptoLastInj_m(0),
      globalPartPerNode_m(nullptr),
      dist_m(nullptr),
      dcBeam_m(false),
      periodLength_m(Physics::c / 1e9),
      pbase_m(pb)
{
    setup(pb);
}

/*
 * Bunch common member functions
 */

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::getIfBeamEmitting() {
    if (dist_m != nullptr) {
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
    for (int i = 0; i < numberOfEnergyBins; i++)
        binemitted_m[i] = 0;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveEnergyBins() {

    if (dist_m != nullptr)
        return dist_m->getNumberOfEnergyBins() > 0;
    else
        return false;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::switchToUnitlessPositions(bool use_dt_per_particle) {

    if (unit_state_ == unitless)
        throw SwitcherError("PartBunch::switchToUnitlessPositions",
                            "Cannot make a unitless PartBunch unitless");

    bool hasToReset = false;
    if (!R.isDirty()) hasToReset = true;

    for (size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt[i];

        R[i] /= Vector_t(Physics::c * dt);
    }

    unit_state_ = unitless;

    if (hasToReset) R.resetDirtyFlag();
}

//FIXME: unify methods, use convention that all particles have own dt
template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::switchOffUnitlessPositions(bool use_dt_per_particle) {

    if (unit_state_ == units)
        throw SwitcherError("PartBunch::switchOffUnitlessPositions",
                            "Cannot apply units twice to PartBunch");

    bool hasToReset = false;
    if (!R.isDirty()) hasToReset = true;

    for (size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt[i];

        R[i] *= Vector_t(Physics::c * dt);
    }

    unit_state_ = units;

    if (hasToReset) R.resetDirtyFlag();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::do_binaryRepart() {
    // do nothing here
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setDistribution(Distribution* d,
                                            std::vector<Distribution*> addedDistributions,
                                            size_t& np) {
    Inform m("setDistribution " );
    dist_m = d;
    dist_m->createOpalT(this, addedDistributions, np);

//    if (Options::cZero)
//        dist_m->create(this, addedDistributions, np / 2);
//    else
//        dist_m->create(this, addedDistributions, np);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setDistribution(Distribution* d,
                                            size_t numberOfParticles,
                                            double current,
                                            const Beamline& bl) {
    dist_m = d;
    dist_m->createOpalCycl(this, numberOfParticles, current, bl);
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::isGridFixed() const {
    return fixed_grid;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::hasBinning() const {
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
    if (dist_m != nullptr)
        return dist_m->getIfDistEmitting();
    else
        return false;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveBins() const {
    if (pbin_m != nullptr)
        return pbin_m->weHaveBins();
    else
        return false;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBins* pbin) {
    pbin_m = pbin;
    *gmsg << *pbin_m << endl;
    setEnergyBins(pbin_m->getNBins());
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBinsCyc* pbin) {
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
    pbin_m = nullptr;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getLastemittedBin() {
    if (pbin_m != nullptr)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalBinCount(size_t num, int bin) {
    if (pbin_m != nullptr) {
        pbin_m->setPartNum(bin, num);
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas() {

    const int emittedBins = dist_m->getNumberOfEnergyBins();

    size_t s = 0;

    for (int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;

    for (unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += std::sqrt(1.0 + dot(this->P[n], this->P[n]));

    std::unique_ptr<size_t[]> particlesInBin(new size_t[emittedBins]);
    reduce(bingamma_m.get(), bingamma_m.get() + emittedBins, bingamma_m.get(), OpAddAssign());
    reduce(binemitted_m.get(), binemitted_m.get() + emittedBins, particlesInBin.get(), OpAddAssign());
    for (int i = 0; i < emittedBins; i++) {
        size_t &pInBin = particlesInBin[i];
        if (pInBin != 0) {
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


    if (s != getTotalNum() && !OpalData::getInstance()->hasGlobalGeometry())
        ERRORMSG("sum(Bins)= " << s << " != sum(R)= " << getTotalNum() << endl;);

    if (emittedBins >= 2) {
        for (int i = 1; i < emittedBins; i++) {
            if (binemitted_m[i - 1] != 0 && binemitted_m[i] != 0)
                INFOMSG(level2 << "d(gamma)= " << 100.0 * std::abs(bingamma_m[i - 1] - bingamma_m[i]) / bingamma_m[i] << " [%] "
                        << "between bin " << i - 1 << " and " << i << endl);
        }
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas_cycl() {

    const int emittedBins = pbin_m->getLastemittedBin();

    for (int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;
    for (unsigned int n = 0; n < getLocalNum(); n++) {
        if ( this->Bin[n] > -1 ) {
            bingamma_m[this->Bin[n]] += std::sqrt(1.0 + dot(this->P[n], this->P[n]));
        }
    }

    allreduce(*bingamma_m.get(), emittedBins, std::plus<double>());

    for (int i = 0; i < emittedBins; i++) {
        if (pbin_m->getTotalNumPerBin(i) > 0) {
            bingamma_m[i] /= pbin_m->getTotalNumPerBin(i);
        } else {
            bingamma_m[i] = 0.0;
        }
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

    for (unsigned long k = 0; k < getLocalNum(); ++ k)
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
                                            std::vector<double>& lineDensity,
                                            std::pair<double, double>& meshInfo) {
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
    //if (!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;

    if (!isGridFixed()) {
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
        for (int i = 0; i < dimIdx; i++) {
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

        if (getIfBeamEmitting() && dist_m != nullptr) {
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

        if (hr_m[0] * hr_m[1] * hr_m[2] <= 0) {
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
void PartBunchBase<T, Dim>::boundp_destroyCycl() {

    Vector_t len;
    const int dimIdx = 3;
    IpplTimings::startTimer(boundpTimer_m);

    std::unique_ptr<size_t[]> countLost;
    if (weHaveBins()) {
        const int tempN = pbin_m->getLastemittedBin();
        countLost = std::unique_ptr<size_t[]>(new size_t[tempN]);
        for (int ii = 0; ii < tempN; ii++) countLost[ii] = 0;
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
        Vector_t rmean = momentsComputer_m.getMeanPosition();
        Vector_t rrms = momentsComputer_m.getStandardDeviationPosition();
        if(checkfactor < 0) {
            checkfactor *= -1;
            if (len[0] > checkfactor * rrms[0] ||
                len[1] > checkfactor * rrms[1] ||
                len[2] > checkfactor * rrms[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the distance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (std::abs(R[ii](0) - rmean(0)) > checkfactor * rrms[0] ||
                        std::abs(R[ii](1) - rmean(1)) > checkfactor * rrms[1] ||
                        std::abs(R[ii](2) - rmean(2)) > checkfactor * rrms[2])
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
            if (len[0] > checkfactor * rrms[0] ||
                len[2] > checkfactor * rrms[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the distance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (std::abs(R[ii](0) - rmean(0)) > checkfactor * rrms[0] ||
                        std::abs(R[ii](2) - rmean(2)) > checkfactor * rrms[2])
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

    for (int i = 0; i < dimIdx; i++) {
        double length = std::abs(rmax_m[i] - rmin_m[i]);
        rmax_m[i] += dh_m * length;
        rmin_m[i] -= dh_m * length;
        hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
    }

    // rescale mesh
    this->updateFields(hr_m, rmin_m);

    if (weHaveBins()) {
        pbin_m->updatePartInBin_cyc(countLost.get());
    }

    /* we also need to update the number of particles per bunch
     * expensive since does an allreduce!
     */
    countTotalNumPerBunch();

    IpplTimings::startTimer(boundpUpdateTimer_m);
    update();
    IpplTimings::stopTimer(boundpUpdateTimer_m);

    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::boundp_destroyT() {

    this->updateDomainLength(nr_m);

    std::vector<size_t> tmpbinemitted;

    boundp();

    size_t ne = 0;
    const size_t localNum = getLocalNum();

    double rzmean = momentsComputer_m.getMeanPosition()(2);
    double rzrms = momentsComputer_m.getStandardDeviationPosition()(2);
    const bool haveEnergyBins = weHaveEnergyBins();
    if (haveEnergyBins) {
        tmpbinemitted.resize(getNumberOfEnergyBins(), 0.0);
    }
    for (unsigned int i = 0; i < localNum; i++) {
        if (Bin[i] < 0 || (Options::remotePartDel > 0 && std::abs(R[i](2) - rzmean) < Options::remotePartDel * rzrms)) {
            ne++;
            destroy(1, i);
        } else if (haveEnergyBins) {
            tmpbinemitted[Bin[i]]++;
        }
    }

    boundp();

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if (haveEnergyBins) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for (int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    reduce(ne, ne, OpAddAssign());
    return ne;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::destroyT() {

    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = getLocalNum();
    const size_t totalNum = getTotalNum();
    size_t ne = 0;

    if (weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for (int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for (unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                destroy(1, i);
                ++ ne;
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        for (size_t i = 0; i < localNum; i++) {
            if ((Bin[i] < 0)) {
                ne++;
                destroy(1, i);
            }
        }
    }

    if (ne > 0) {
        performDestroy(true);
    }

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if (weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for (int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    size_t newTotalNum = getLocalNum();
    reduce(newTotalNum, newTotalNum, OpAddAssign());

    setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPx(int /*i*/) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPy(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPz(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPx0(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPy0(int) {
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
double PartBunchBase<T, Dim>::getX0(int /*i*/) {
    return 0.0;
}

//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getY0(int /*i*/) {
    return 0.0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setZ(int /*i*/, double /*zcoo*/) {
};


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::get_bounds(Vector_t& rmin, Vector_t& rmax) const {

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
void PartBunchBase<T, Dim>::getLocalBounds(Vector_t& rmin, Vector_t& rmax) const {
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
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
std::pair<Vector_t, double> PartBunchBase<T, Dim>::getLocalBoundingSphere() {
    Vector_t rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::push_back(OpalParticle const& particle) {
    Inform msg("PartBunch ");

    size_t i = getLocalNum();
    create(1);

    R[i] = particle.getR();
    P[i] = particle.getP();

    update();
    msg << "Created one particle i= " << i << endl;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setParticle(FVector<double, 6> z, int ii) {
    R[ii](0) = z[0];
    P[ii](0) = z[1];
    R[ii](1) = z[2];
    P[ii](1) = z[3];
    R[ii](2) = z[4];
    P[ii](2) = z[5];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setParticle(OpalParticle const& particle, int ii) {
    R[ii] = particle.getR();
    P[ii] = particle.getP();
}


template <class T, unsigned Dim>
OpalParticle PartBunchBase<T, Dim>::getParticle(int ii) {
    OpalParticle particle(ID[ii],
                          Vector_t(R[ii](0), R[ii](1), 0), P[ii],
                          R[ii](2), Q[ii], M[ii]);
    return particle;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::maximumAmplitudes(const FMatrix<double, 6, 6>& D,
                                              double& axmax, double& aymax) {
    axmax = aymax = 0.0;
    OpalParticle particle;

    for (unsigned int ii = 0; ii < getLocalNum(); ii++) {

        particle = getParticle(ii);
        FVector<double, 6> vec({particle.getX(), particle.getPx(),
                                particle.getY(), particle.getPy(),
                                particle.getZ(), particle.getPz()});
        FVector<double, 6> result;
        result = D * vec;
        // double xnor =
        //     D(0, 0) * part.getX()  + D(0, 1) * part.getPx() + D(0, 2) * part.getY() +
        //     D(0, 3) * part.getPy() + D(0, 4) * part.getL()  + D(0, 5) * part.getPLon();
        // double pxnor =
        //     D(1, 0) * part.getX()  + D(1, 1) * part.getPx() + D(1, 2) * part.getY() +
        //     D(1, 3) * part.getPy() + D(1, 4) * part.getL()  + D(1, 5) * part.getPLon();
        // double ynor =
        //     D(2, 0) * part.getX()  + D(2, 1) * part.getPx() + D(2, 2) * part.getY() +
        //     D(2, 3) * part.getPy() + D(2, 4) * part.getL()  + D(2, 5) * part.getPLon();
        // double pynor =
        //     D(3, 0) * part.getX()  + D(3, 1) * part.getPx() + D(3, 2) * part.getY() +
        //     D(3, 3) * part.getPy() + D(3, 4) * part.getL()  + D(3, 5) * part.getPLon();

        axmax = std::max(axmax, (std::pow(result[0], 2) + std::pow(result[1], 2)));
        aymax = std::max(aymax, (std::pow(result[2], 2) + std::pow(result[3], 2)));
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
void PartBunchBase<T, Dim>::incrementT() {
    t_m += dt_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getT() const {
    return t_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_sPos() const {
    return spos_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_sPos(double s) {
    spos_m = s;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_gamma() const {
    return momentsComputer_m.getMeanGamma();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_meanKineticEnergy() const {
    return momentsComputer_m.getMeanKineticEnergy();
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
    return momentsComputer_m.getMeanPosition();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rrms() const {
    return momentsComputer_m.getStandardDeviationPosition();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rprms() const {
    return momentsComputer_m.getStandardDeviationRP();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_rmean() const {
    return momentsComputer_m.getMeanPosition();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_prms() const {
    return momentsComputer_m.getStandardDeviationMomentum();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_pmean() const {
    return momentsComputer_m.getMeanMomentum();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_pmean_Distribution() const {
    if (dist_m)// && dist_m->getType() != DistrTypeT::FROMFILE)
        return dist_m->get_pmean();

    double gamma = 0.1 / getM() + 1; // set default 0.1 eV
    return Vector_t(0, 0, std::sqrt(std::pow(gamma, 2) - 1));
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_emit() const {
    return momentsComputer_m.getGeometricEmittance();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_norm_emit() const {
    return momentsComputer_m.getNormalizedEmittance();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_halo() const {
    return momentsComputer_m.getHalo();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_68Percentile() const {
    return momentsComputer_m.get68Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_95Percentile() const {
    return momentsComputer_m.get95Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_99Percentile() const {
    return momentsComputer_m.get99Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_99_99Percentile() const {
    return momentsComputer_m.get99_99Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_normalizedEps_68Percentile() const {
    return momentsComputer_m.getNormalizedEmittance68Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_normalizedEps_95Percentile() const {
    return momentsComputer_m.getNormalizedEmittance95Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_normalizedEps_99Percentile() const {
    return momentsComputer_m.getNormalizedEmittance99Percentile();
}

template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_normalizedEps_99_99Percentile() const {
    return momentsComputer_m.getNormalizedEmittance99_99Percentile();
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_Dx() const {
    return momentsComputer_m.getDx();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_Dy() const {
    return momentsComputer_m.getDy();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_DDx() const {
    return momentsComputer_m.getDDx();
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::get_DDy() const {
    return momentsComputer_m.getDDy();
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_hr() const {
    return hr_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::set_meshEnlargement(double dh) {
    dh_m = dh;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::gatherLoadBalanceStatistics() {

    for (int i = 0; i < Ippl::getNodes(); i++)
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

    IpplTimings::startTimer(statParamTimer_m);
    get_bounds(rmin_m, rmax_m);
    momentsComputer_m.compute(*this);
    IpplTimings::stopTimer(statParamTimer_m);
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
    if (getTotalNum() != 0)
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
    massPerParticle_m = mass;
    if (getTotalNum() != 0)
        M = mass;
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setMassZeroPart(double mass) {
    massPerParticle_m = mass;
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
double PartBunchBase<T, Dim>::getMassPerParticle() const {
    return massPerParticle_m;
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setSolver(FieldSolver *fs) {
    fs_m = fs;
    fs_m->initSolver(this);

    /**
       CAN not re-inizialize ParticleLayout
       this is an IPPL issue
     */
    if (!OpalData::getInstance()->hasBunchAllocated()) {
        this->initialize(fs_m->getFieldLayout());
//         this->setMesh(fs_m->getMesh());
//         this->setFieldLayout(fs_m->getFieldLayout());
    }
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::hasFieldSolver() {
    if (fs_m)
        return fs_m->hasValidSolver();
    else
        return false;
}


/// \brief Return the fieldsolver type if we have a fieldsolver
template <class T, unsigned Dim>
FieldSolverType PartBunchBase<T, Dim>::getFieldSolverType() const {
    if (fs_m) {
        return fs_m->getFieldSolverType();
    } else {
        return FieldSolverType::NONE;
    }
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
void PartBunchBase<T, Dim>::setNumBunch(short n) {
    numBunch_m = n;
    bunchTotalNum_m.resize(n);
    bunchLocalNum_m.resize(n);
}


template <class T, unsigned Dim>
short PartBunchBase<T, Dim>::getNumBunch() const {
    return numBunch_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setTotalNumPerBunch(size_t totalnum, short n) {
    PAssert_LT(n, (short)bunchTotalNum_m.size());
    bunchTotalNum_m[n] = totalnum;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getTotalNumPerBunch(short n) const {
    PAssert_LT(n, (short)bunchTotalNum_m.size());
    return bunchTotalNum_m[n];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalNumPerBunch(size_t localnum, short n) {
    PAssert_LT(n, (short)bunchLocalNum_m.size());
    bunchLocalNum_m[n] = localnum;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getLocalNumPerBunch(short n) const {
    PAssert_LT(n, (short)bunchLocalNum_m.size());
    return bunchLocalNum_m[n];
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::countTotalNumPerBunch() {
    bunchTotalNum_m.clear();
    bunchTotalNum_m.resize(numBunch_m);
    bunchLocalNum_m.clear();
    bunchLocalNum_m.resize(numBunch_m);

    for (size_t i = 0; i < this->getLocalNum(); ++i) {
        PAssert_LT(this->bunchNum[i], numBunch_m);
        ++bunchLocalNum_m[this->bunchNum[i]];
    }

    allreduce(bunchLocalNum_m.data(), bunchTotalNum_m.data(),
              bunchLocalNum_m.size(), std::plus<size_t>());

    size_t totalnum = std::accumulate(bunchTotalNum_m.begin(),
                                      bunchTotalNum_m.end(), 0);

    if ( totalnum != this->getTotalNum() )
        throw OpalException("PartBunchBase::countTotalNumPerBunch()",
                            "Sum of total number of particles per bunch (" +
                            std::to_string(totalnum) + ") != total number of particles (" +
                            std::to_string(this->getTotalNum()) + ")");
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
int PartBunchBase<T, Dim>::getSteptoLastInj() const {
    return SteptoLastInj_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calcMeanPhi() {

    const int emittedBins = pbin_m->getLastemittedBin();
    double phi[emittedBins];
    double px[emittedBins];
    double py[emittedBins];
    double meanPhi = 0.0;

    for (int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = 0.0;
        px[ii] = 0.0;
        py[ii] = 0.0;
    }

    for (unsigned int ii = 0; ii < getLocalNum(); ii++) {
        px[Bin[ii]] += P[ii](0);
        py[Bin[ii]] += P[ii](1);
    }

    reduce(px, px + emittedBins, px, OpAddAssign());
    reduce(py, py + emittedBins, py, OpAddAssign());
    for (int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = calculateAngle(px[ii], py[ii]);
        meanPhi += phi[ii];
        INFOMSG("Bin " << ii  << "  mean phi = " << phi[ii] * Units::rad2deg - 90.0 << "[degree]" << endl);
    }

    meanPhi /= emittedBins;

    INFOMSG("mean phi of all particles " <<  meanPhi * Units::rad2deg - 90.0 << "[degree]" << endl);

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
    for (int ii = 0; ii < maxbin; ii++) partInBin[ii] = 0;

    double pMin0 = 1.0e9;
    double pMin = 0.0;
    double maxbinIndex = 0;

    for (unsigned long int n = 0; n < getLocalNum(); n++) {
        double temp_betagamma = std::sqrt(std::pow(P[n](0), 2) + std::pow(P[n](1), 2));
        if (pMin0 > temp_betagamma)
            pMin0 = temp_betagamma;
    }
    reduce(pMin0, pMin, OpMinAssign());
    INFOMSG("minimal beta*gamma = " << pMin << endl);

    double asinh0 = std::asinh(pMin);
    for (unsigned long int n = 0; n < getLocalNum(); n++) {

        double temp_betagamma = std::sqrt(std::pow(P[n](0), 2) + std::pow(P[n](1), 2));
        int itsBinID = std::floor((std::asinh(temp_betagamma) - asinh0) / eta + 1.0E-6);
        Bin[n] = itsBinID;
        if (maxbinIndex < itsBinID) {
            maxbinIndex = itsBinID;
        }

        if (itsBinID >= maxbin) {
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
ParticleOrigin PartBunchBase<T, Dim>::getPOrigin() const {
    return refPOrigin_m;
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPOrigin(ParticleOrigin origin) {
    refPOrigin_m = origin;
}


template <class T, unsigned Dim>
ParticleType PartBunchBase<T, Dim>::getPType() const {
    return refPType_m;
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPType(const std::string& type) {
    refPType_m =  ParticleProperties::getParticleType(type);
}


template <class T, unsigned Dim>
DistributionType PartBunchBase<T, Dim>::getDistType() const {
    return dist_m->getType();
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
double PartBunchBase<T, Dim>::getdE() const {
    return momentsComputer_m.getStdKineticEnergy();
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
double PartBunchBase<T, Dim>::getGamma(int /*i*/) {
    return 0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getBeta(int /*i*/) {
    return 0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::actT() {
    // do nothing here
};


template <class T, unsigned Dim>
const PartData* PartBunchBase<T, Dim>::getReference() const {
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
    momentsComputer_m.computeMeanKineticEnergy(*this);
}

template <class T, unsigned Dim>
Inform& PartBunchBase<T, Dim>::print(Inform& os) {

    if (getTotalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();

        double pathLength = get_sPos();

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H ********************************************************* \n";
        os << "* NP              = " << getTotalNum() << "\n";
        os << "* Qtot            = " << std::setw(17) << Util::getChargeString(std::abs(sum(Q))) << "         "
        << "Qi    = "             << std::setw(17) << Util::getChargeString(std::abs(qi_m)) << "\n";
        os << "* Ekin            = " << std::setw(17) << Util::getEnergyString(get_meanKineticEnergy()) << "         "
           << "dEkin = "             << std::setw(17) << Util::getEnergyString(getdE()) << "\n";
        os << "* rmax            = " << Util::getLengthString(rmax_m, 5) << "\n";
        os << "* rmin            = " << Util::getLengthString(rmin_m, 5) << "\n";
        if (getTotalNum() >= 2) { // to suppress Nans
            os << "* rms beam size   = " << Util::getLengthString(get_rrms(), 5) << "\n";
            os << "* rms momenta     = " << std::setw(12) << std::setprecision(5) << get_prms() << " [beta gamma]\n";
            os << "* mean position   = " << Util::getLengthString(get_rmean(), 5) << "\n";
            os << "* mean momenta    = " << std::setw(12) << std::setprecision(5) << get_pmean() << " [beta gamma]\n";
            os << "* rms emittance   = " << std::setw(12) << std::setprecision(5) << get_emit() << " (not normalized)\n";
            os << "* rms correlation = " << std::setw(12) << std::setprecision(5) << get_rprms() << "\n";
        }
        os << "* hr              = " << Util::getLengthString(get_hr(), 5) << "\n";
        os << "* dh              = " << std::setw(13) << std::setprecision(5) << dh_m * 100 << " [%]\n";
        os << "* t               = " << std::setw(17) << Util::getTimeString(getT()) << "         "
           << "dT    = "             << std::setw(17) << Util::getTimeString(getdT()) << "\n";
        os << "* spos            = " << std::setw(17) << Util::getLengthString(pathLength) << "\n";
        os << "* ********************************************************************************** " << endl;
        os.flags(ff);
    }
    return os;
}

// angle range [0~2PI) degree
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    return thetaXY >= 0 ? thetaXY : thetaXY + Physics::two_pi;
}


template <class T, unsigned Dim>
Inform& operator<<(Inform &os, PartBunchBase<T, Dim>& p) {
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
void PartBunchBase<T, Dim>::resetInterpolationCache(bool /*clearCache*/) {

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
    std::swap(POrigin[i], POrigin[j]);
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
void PartBunchBase<T, Dim>::updateFields(const Vector_t& /*hr*/, const Vector_t& /*origin*/) {
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
    pb->addAttribute(POrigin);
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

    pmsg_m.release();
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getTotalNum() const {
    return pbase_m->getTotalNum();
}

template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getLocalNum() const {
    return pbase_m->getLocalNum();
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getDestroyNum() const {
    return pbase_m->getDestroyNum();
}

template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getGhostNum() const {
    return pbase_m->getGhostNum();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setTotalNum(size_t n) {
    pbase_m->setTotalNum(n);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setLocalNum(size_t n) {
    pbase_m->setLocalNum(n);
}

template <class T, unsigned Dim>
ParticleLayout<T, Dim> & PartBunchBase<T, Dim>::getLayout() {
    return pbase_m->getLayout();
}

template <class T, unsigned Dim>
const ParticleLayout<T, Dim>& PartBunchBase<T, Dim>::getLayout() const {
    return pbase_m->getLayout();
}

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::getUpdateFlag(UpdateFlags_t f) const {
    return pbase_m->getUpdateFlag(f);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setUpdateFlag(UpdateFlags_t f, bool val) {
    pbase_m->setUpdateFlag(f, val);
}

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::singleInitNode() const {
    return pbase_m->singleInitNode();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetID() {
    pbase_m->resetID();
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::update() {
    try {
        pbase_m->update();
    } catch (const IpplException& ex) {
        throw OpalException(ex.where(), ex.what());
    }
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::update(const ParticleAttrib<char>& canSwap) {
    try {
        pbase_m->update(canSwap);
    } catch (const IpplException& ex) {
        throw OpalException(ex.where(), ex.what());
    }
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::createWithID(unsigned id) {
    pbase_m->createWithID(id);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::create(size_t M) {
    pbase_m->create(M);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::globalCreate(size_t np) {
    pbase_m->globalCreate(np);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::destroy(size_t M, size_t I, bool doNow) {
    pbase_m->destroy(M, I, doNow);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::performDestroy(bool updateLocalNum) {
    pbase_m->performDestroy(updateLocalNum);
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::ghostDestroy(size_t M, size_t I) {
    pbase_m->ghostDestroy(M, I);
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
        rpmean(2*i)= get_rmean()(i);
        rpmean((2*i)+1)= get_pmean()(i);
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
