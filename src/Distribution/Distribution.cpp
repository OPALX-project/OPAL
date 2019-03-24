// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//   The class for the OPAL Distribution command.
//
// ------------------------------------------------------------------------

#include "Distribution/Distribution.h"
#include "Distribution/SigmaGenerator.h"
#include "AbsBeamline/SpecificElementVisitor.h"

#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <limits>

#include "AbstractObjects/Expressions.h"
#include "Utilities/Options.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBunchBase.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Algorithms/PartBinsCyc.h"
#include "BasicActions/Option.h"
#include "Distribution/LaserProfile.h"
#include "Elements/OpalBeamline.h"
#include "AbstractObjects/BeamSequence.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPC.h"
#include "Utilities/Util.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <sys/time.h>

#include <boost/regex.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

extern Inform *gmsg;

constexpr double SMALLESTCUTOFF = 1e-12;

namespace {
    SymTenzor<double, 6> getUnit6x6() {
        SymTenzor<double, 6> unit6x6;
        for (unsigned int i = 0; i < 6u; ++ i) {
            unit6x6(i,i) = 1.0;
        }
        return unit6x6;
    }
}

//
// Class Distribution
// ------------------------------------------------------------------------

Distribution::Distribution():
    Definition( Attrib::Legacy::Distribution::SIZE, "DISTRIBUTION",
                "The DISTRIBUTION statement defines data for the 6D particle distribution."),
    distrTypeT_m(DistrTypeT::NODIST),
    numberOfDistributions_m(1),
    emitting_m(false),
    emissionModel_m(EmissionModelT::NONE),
    tEmission_m(0.0),
    tBin_m(0.0),
    currentEmissionTime_m(0.0),
    currentEnergyBin_m(1),
    currentSampleBin_m(0),
    numberOfEnergyBins_m(0),
    numberOfSampleBins_m(0),
    energyBins_m(NULL),
    energyBinHist_m(NULL),
    randGen_m(NULL),
    pTotThermal_m(0.0),
    pmean_m(0.0),
    cathodeWorkFunc_m(0.0),
    laserEnergy_m(0.0),
    cathodeFermiEnergy_m(0.0),
    cathodeTemp_m(0.0),
    emitEnergyUpperLimit_m(0.0),
    totalNumberParticles_m(0),
    totalNumberEmittedParticles_m(0),
    avrgpz_m(0.0),
    inputMoUnits_m(InputMomentumUnitsT::NONE),
    sigmaTRise_m(0.0),
    sigmaTFall_m(0.0),
    tPulseLengthFWHM_m(0.0),
    correlationMatrix_m(getUnit6x6()),
    laserProfileFileName_m(""),
    laserImageName_m(""),
    laserIntensityCut_m(0.0),
    laserProfile_m(NULL),
    darkCurrentParts_m(0),
    darkInwardMargin_m(0.0),
    eInitThreshold_m(0.0),
    workFunction_m(0.0),
    fieldEnhancement_m(0.0),
    fieldThrFN_m(0.0),
    maxFN_m(0),
    paraFNA_m(0.0),
    paraFNB_m(0.0),
    paraFNY_m(0.0),
    paraFNVYSe_m(0.0),
    paraFNVYZe_m(0.0),
    secondaryFlag_m(0),
    ppVw_m(0.0),
    vVThermal_m(0.0),
    I_m(0.0),
    E_m(0.0)
{
    setAttributes();

    Distribution *defaultDistribution = clone("UNNAMED_Distribution");
    defaultDistribution->builtin = true;

    try {
        OpalData::getInstance()->define(defaultDistribution);
    } catch(...) {
        delete defaultDistribution;
    }

    // setFieldEmissionParameters();

    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);
}
/**
 *
 *
 * @param name
 * @param parent
 */
Distribution::Distribution(const std::string &name, Distribution *parent):
    Definition(name, parent),
    distT_m(parent->distT_m),
    distrTypeT_m(DistrTypeT::NODIST),
    numberOfDistributions_m(parent->numberOfDistributions_m),
    emitting_m(parent->emitting_m),
    particleRefData_m(parent->particleRefData_m),
    addedDistributions_m(parent->addedDistributions_m),
    particlesPerDist_m(parent->particlesPerDist_m),
    emissionModel_m(parent->emissionModel_m),
    tEmission_m(parent->tEmission_m),
    tBin_m(parent->tBin_m),
    currentEmissionTime_m(parent->currentEmissionTime_m),
    currentEnergyBin_m(parent->currentEmissionTime_m),
    currentSampleBin_m(parent->currentSampleBin_m),
    numberOfEnergyBins_m(parent->numberOfEnergyBins_m),
    numberOfSampleBins_m(parent->numberOfSampleBins_m),
    energyBins_m(NULL),
    energyBinHist_m(NULL),
    randGen_m(NULL),
    pTotThermal_m(parent->pTotThermal_m),
    pmean_m(parent->pmean_m),
    cathodeWorkFunc_m(parent->cathodeWorkFunc_m),
    laserEnergy_m(parent->laserEnergy_m),
    cathodeFermiEnergy_m(parent->cathodeFermiEnergy_m),
    cathodeTemp_m(parent->cathodeTemp_m),
    emitEnergyUpperLimit_m(parent->emitEnergyUpperLimit_m),
    totalNumberParticles_m(parent->totalNumberParticles_m),
    totalNumberEmittedParticles_m(parent->totalNumberEmittedParticles_m),
    xDist_m(parent->xDist_m),
    pxDist_m(parent->pxDist_m),
    yDist_m(parent->yDist_m),
    pyDist_m(parent->pyDist_m),
    tOrZDist_m(parent->tOrZDist_m),
    pzDist_m(parent->pzDist_m),
    xWrite_m(parent->xWrite_m),
    pxWrite_m(parent->pxWrite_m),
    yWrite_m(parent->yWrite_m),
    pyWrite_m(parent->pyWrite_m),
    tOrZWrite_m(parent->tOrZWrite_m),
    pzWrite_m(parent->pzWrite_m),
    avrgpz_m(parent->avrgpz_m),
    inputMoUnits_m(parent->inputMoUnits_m),
    sigmaTRise_m(parent->sigmaTRise_m),
    sigmaTFall_m(parent->sigmaTFall_m),
    tPulseLengthFWHM_m(parent->tPulseLengthFWHM_m),
    sigmaR_m(parent->sigmaR_m),
    sigmaP_m(parent->sigmaP_m),
    cutoffR_m(parent->cutoffR_m),
    cutoffP_m(parent->cutoffP_m),
    correlationMatrix_m(parent->correlationMatrix_m),
    laserProfileFileName_m(parent->laserProfileFileName_m),
    laserImageName_m(parent->laserImageName_m),
    laserIntensityCut_m(parent->laserIntensityCut_m),
    laserProfile_m(NULL),
    darkCurrentParts_m(parent->darkCurrentParts_m),
    darkInwardMargin_m(parent->darkInwardMargin_m),
    eInitThreshold_m(parent->eInitThreshold_m),
    workFunction_m(parent->workFunction_m),
    fieldEnhancement_m(parent->fieldEnhancement_m),
    fieldThrFN_m(parent->fieldThrFN_m),
    maxFN_m(parent->maxFN_m),
    paraFNA_m(parent-> paraFNA_m),
    paraFNB_m(parent-> paraFNB_m),
    paraFNY_m(parent-> paraFNY_m),
    paraFNVYSe_m(parent-> paraFNVYSe_m),
    paraFNVYZe_m(parent-> paraFNVYZe_m),
    secondaryFlag_m(parent->secondaryFlag_m),
    ppVw_m(parent->ppVw_m),
    vVThermal_m(parent->vVThermal_m),
    I_m(parent->I_m),
    E_m(parent->E_m),
    tRise_m(parent->tRise_m),
    tFall_m(parent->tFall_m),
    sigmaRise_m(parent->sigmaRise_m),
    sigmaFall_m(parent->sigmaFall_m),
    cutoff_m(parent->cutoff_m)
{
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);
}

Distribution::~Distribution() {

    delete energyBins_m;
    gsl_histogram_free(energyBinHist_m);
    gsl_rng_free(randGen_m);
    delete laserProfile_m;
}

/*
  void Distribution::printSigma(SigmaGenerator<double,unsigned int>::matrix_type& M, Inform& out) {
  for (int i=0; i<M.size1(); ++i) {
  for (int j=0; j<M.size2(); ++j) {
  *gmsg  << M(i,j) << " ";
  }
  *gmsg << endl;
  }
  }
*/

/**
 * Calculate the local number of particles evenly and adjust node 0
 * such that n is matched exactly.
 * @param n total number of particles
 * @return n / #cores
 * @param
 */

size_t Distribution::getNumOfLocalParticlesToCreate(size_t n) {

    size_t locNumber = n / Ippl::getNodes();

    // make sure the total number is exact
    size_t remainder  = n % Ippl::getNodes();
    size_t myNode = Ippl::myNode();
    if (myNode < remainder)
        ++ locNumber;

    return locNumber;
}

/// Distribution can only be replaced by another distribution.
bool Distribution::canReplaceBy(Object *object) {
    return dynamic_cast<Distribution *>(object) != 0;
}

Distribution *Distribution::clone(const std::string &name) {
    return new Distribution(name, this);
}

void Distribution::execute() {
}

void Distribution::update() {
}

void Distribution::create(size_t &numberOfParticles, double massIneV) {

    /*
     * If Options::cZero is true than we reflect generated distribution
     * to ensure that the transverse averages are 0.0.
     *
     * For now we just cut the number of generated particles in half.
     */
    size_t numberOfLocalParticles = numberOfParticles;
    if (Options::cZero && distrTypeT_m != DistrTypeT::FROMFILE) {
        numberOfLocalParticles = (numberOfParticles + 1) / 2;
    }

    size_t mySeed = Options::seed;

    if (Options::seed == -1) {
        struct timeval tv;
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec + tv.tv_usec;
    }

    if (Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE])) {
        numberOfLocalParticles = getNumOfLocalParticlesToCreate(numberOfLocalParticles);
        *gmsg << level2 << "* Generation of distribution with seed = " << mySeed << " + core_id\n"
              << "* is scalable with number of particles and cores." << endl;
        mySeed += Ippl::myNode();
    } else {
        *gmsg << level2 << "* Generation of distribution with seed = " << mySeed << "\n"
              << "* isn't scalable with number of particles and cores." << endl;
    }

    gsl_rng_set(randGen_m, mySeed);

    setFieldEmissionParameters();

    switch (distrTypeT_m) {

    case DistrTypeT::MATCHEDGAUSS:
        createMatchedGaussDistribution(numberOfLocalParticles, massIneV);
        break;
    case DistrTypeT::FROMFILE:
        createDistributionFromFile(numberOfParticles, massIneV);
        break;
    case DistrTypeT::GAUSS:
        createDistributionGauss(numberOfLocalParticles, massIneV);
        break;
    case DistrTypeT::BINOMIAL:
        createDistributionBinomial(numberOfLocalParticles, massIneV);
        break;
    case DistrTypeT::FLATTOP:
    case DistrTypeT::GUNGAUSSFLATTOPTH:
    case DistrTypeT::ASTRAFLATTOPTH:
        createDistributionFlattop(numberOfLocalParticles, massIneV);
        break;
    default:
        INFOMSG("Distribution unknown." << endl;);
        break;
    }

    if (emitting_m) {

        unsigned int numAdditionalRNsPerParticle = 0;
        if (emissionModel_m == EmissionModelT::ASTRA ||
            distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH ||
            distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH) {

            numAdditionalRNsPerParticle = 2;
        } else if (emissionModel_m == EmissionModelT::NONEQUIL) {
            if (Options::cZero) {
                numAdditionalRNsPerParticle = 40;
            } else {
                numAdditionalRNsPerParticle = 20;
            }
        }

        int saveProcessor = -1;
        const int myNode = Ippl::myNode();
        const int numNodes = Ippl::getNodes();
        const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

        for (size_t partIndex = 0; partIndex < numberOfLocalParticles; ++ partIndex) {

            // Save to each processor in turn.
            ++ saveProcessor;
            if (saveProcessor >= numNodes)
                saveProcessor = 0;

            if (scalable || myNode == saveProcessor) {
                std::vector<double> rns;
                for (unsigned int i = 0; i < numAdditionalRNsPerParticle; ++ i) {
                    double x = gsl_rng_uniform(randGen_m);
                    rns.push_back(x);
                }
                additionalRNs_m.push_back(rns);
            } else {
                for (unsigned int i = 0; i < numAdditionalRNsPerParticle; ++ i) {
                    gsl_rng_uniform(randGen_m);
                }
            }
        }
    }

    // Scale coordinates according to distribution input.
    scaleDistCoordinates();

    // Now reflect particles if Options::cZero is true
    reflectDistribution(numberOfLocalParticles);

    adjustPhaseSpace(massIneV);

    if (Options::seed != -1)
        Options::seed = gsl_rng_uniform_int(randGen_m, gsl_rng_max(randGen_m));

    if (particlesPerDist_m.size() == 0) {
        particlesPerDist_m.push_back(tOrZDist_m.size());
    } else {
        particlesPerDist_m[0] = tOrZDist_m.size();
    }
}

void  Distribution::createPriPart(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg) {

    if (Options::ppdebug) {  // This is Parallel Plate Benchmark.
        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        double vw = this->getVw();
        double vt = this->getvVThermal();
        double f_max = vw / vt * exp(-0.5);
        double test_a = vt / vw;
        double test_asq = test_a * test_a;
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
        if (Ippl::myNode() == 0)
            N_mean += N_extra;
        if (bg.getN() != 0) {
            for (size_t i = 0; i < bg.getN(); i++) {
                if (pc == Ippl::myNode()) {
                    if (count < N_mean) {
                        /*==============Parallel Plate Benchmark=====================================*/
                        double test_s = 1;
                        double f_x = 0;
                        double test_x = 0;
                        while(test_s > f_x) {
                            test_s = gsl_rng_uniform(randGen_m);
                            test_s *= f_max;
                            test_x = gsl_rng_uniform(randGen_m);
                            test_x *= 10 * test_a; //range for normalized emission speed(0,10*test_a);
                            f_x = test_x / test_asq * exp(-test_x * test_x / 2 / test_asq);
                        }
                        double v_emi = test_x * vw;

                        double betaemit = v_emi / Physics::c;
                        double betagamma = betaemit / sqrt(1 - betaemit * betaemit);
                        /*============================================================================ */
                        beam->create(1);
                        if (pc != 0) {
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(Ippl::myNode() * N_mean + count);
                        } else {
                            beam->R[lowMark + count] = bg.getCooridinate(count);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(count);
                        }
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = ParticleType::REGULAR;
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count ++;
                    }
                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
            bg.clearCooridinateArray();
            bg.clearMomentaArray();
            beam->boundp();
        }
        *gmsg << *beam << endl;

    } else {// Normal procedure to create primary particles

        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());

        if (Ippl::myNode() == 0)
            N_mean += N_extra;
        if (bg.getN() != 0) {
            for (size_t i = 0; i < bg.getN(); i++) {

                if (pc == Ippl::myNode()) {
                    if (count < N_mean) {
                        beam->create(1);
                        if (pc != 0)
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra); // node 0 will emit the particle with coordinate ID from 0 to N_mean+N_extra, so other nodes should shift to node_number*N_mean+N_extra
                        else
                            beam->R[lowMark + count] = bg.getCooridinate(count); // for node0 the particle number N_mean =  N_mean + N_extra
                        beam->P[lowMark + count] = Vector_t(0.0);
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = ParticleType::REGULAR;
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count++;

                    }

                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;

            }

        }
        bg.clearCooridinateArray();
        beam->boundp();//fixme if bg.getN()==0?
    }
    *gmsg << *beam << endl;
}

void Distribution::doRestartOpalT(PartBunchBase<double, 3> *beam, size_t Np, int restartStep, H5PartWrapper *dataSource) {

    IpplTimings::startTimer(beam->distrReload_m);

    long numParticles = dataSource->getNumParticles();
    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;
    OpalData::getInstance()->addProblemCharacteristicValue("NP", numParticles);

    numParticles = lastParticle - firstParticle + 1;
    PAssert_GE(numParticles, 0);

    beam->create(numParticles);

    dataSource->readHeader();
    dataSource->readStep(beam, firstParticle, lastParticle);

    beam->boundp();

    double actualT = beam->getT();
    long long ltstep = beam->getLocalTrackStep();
    long long gtstep = beam->getGlobalTrackStep();

    IpplTimings::stopTimer(beam->distrReload_m);

    *gmsg << "Total number of particles in the h5 file= " << beam->getTotalNum() << "\n"
          << "Global step= " << gtstep << "; Local step= " << ltstep << "; "
          << "restart step= " << restartStep << "\n"
          << "time of restart= " << actualT << "; phishift= " << OpalData::getInstance()->getGlobalPhaseShift() << endl;
}

void Distribution::doRestartOpalCycl(PartBunchBase<double, 3> *beam,
                                     size_t Np,
                                     int restartStep,
                                     const int specifiedNumBunch,
                                     H5PartWrapper *dataSource) {

    // h5_int64_t rc;
    IpplTimings::startTimer(beam->distrReload_m);
    INFOMSG("---------------- Start reading hdf5 file----------------" << endl);

    long numParticles = dataSource->getNumParticles();
    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;
    OpalData::getInstance()->addProblemCharacteristicValue("NP", numParticles);

    numParticles = lastParticle - firstParticle + 1;
    PAssert_GE(numParticles, 0);

    beam->create(numParticles);

    dataSource->readHeader();
    dataSource->readStep(beam, firstParticle, lastParticle);

    beam->Q = beam->getChargePerParticle();

    beam->boundp();

    double meanE = static_cast<H5PartWrapperForPC*>(dataSource)->getMeanKineticEnergy();

    const int globalN = beam->getTotalNum();
    INFOMSG("Restart from hdf5 format file " << OpalData::getInstance()->getRestartFileName() << endl);
    INFOMSG("total number of particles = " << globalN << endl);
    INFOMSG("* Restart Energy = " << meanE << " (MeV), Path lenght = " << beam->getLPath() << " (m)" <<  endl);
    INFOMSG("Tracking Step since last bunch injection is " << beam->getSteptoLastInj() << endl);
    INFOMSG(beam->getNumBunch() << " Bunches(bins) exist in this file" << endl);

    double gamma = 1 + meanE / beam->getM() * 1.0e6;
    double beta = sqrt(1.0 - (1.0 / std::pow(gamma, 2.0)));

    INFOMSG("* Gamma = " << gamma << ", Beta = " << beta << endl);

    if (dataSource->predecessorIsSameFlavour()) {
        INFOMSG("Restart from hdf5 file generated by OPAL-cycl" << endl);
        if (specifiedNumBunch > 1) {
            // the allowed maximal bin number is set to 1000
            energyBins_m = new PartBinsCyc(1000, beam->getNumBunch());
            beam->setPBins(energyBins_m);
        }

    } else {
        INFOMSG("Restart from hdf5 file generated by OPAL-t" << endl);

        Vector_t meanR(0.0, 0.0, 0.0);
        Vector_t meanP(0.0, 0.0, 0.0);
        unsigned long int newLocalN = beam->getLocalNum();
        for (unsigned int i = 0; i < newLocalN; ++i) {
            for (int d = 0; d < 3; ++d) {
                meanR(d) += beam->R[i](d);
                meanP(d) += beam->P[i](d);
            }
        }
        reduce(meanR, meanR, OpAddAssign());
        meanR /= Vector_t(globalN);
        reduce(meanP, meanP, OpAddAssign());
        meanP /= Vector_t(globalN);
        INFOMSG("Rmean = " << meanR << "[m], Pmean=" << meanP << endl);

        for (unsigned int i = 0; i < newLocalN; ++i) {
            beam->R[i] -= meanR;
            beam->P[i] -= meanP;
        }
    }

    INFOMSG("---------------Finished reading hdf5 file---------------" << endl);
    IpplTimings::stopTimer(beam->distrReload_m);
}

void Distribution::doRestartOpalE(EnvelopeBunch *beam, size_t Np, int restartStep,
                                  H5PartWrapper *dataSource) {
    IpplTimings::startTimer(beam->distrReload_m);
    int N = dataSource->getNumParticles();
    *gmsg << "total number of slices = " << N << endl;

    beam->distributeSlices(N);
    beam->createBunch();
    long long starti = beam->mySliceStartOffset();
    long long endi = beam->mySliceEndOffset();

    dataSource->readHeader();
    dataSource->readStep(beam, starti, endi);

    beam->setCharge(beam->getChargePerParticle());
    IpplTimings::stopTimer(beam->distrReload_m);
}

Distribution *Distribution::find(const std::string &name) {
    Distribution *dist = dynamic_cast<Distribution *>(OpalData::getInstance()->find(name));

    if (dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}

double Distribution::getTEmission() {
    if (tEmission_m > 0.0) {
        return tEmission_m;
    }

    setDistType();

    tPulseLengthFWHM_m = Attributes::getReal(itsAttr[Attrib::Distribution::TPULSEFWHM]);
    cutoff_m = Attributes::getReal(itsAttr[Attrib::Legacy::Distribution::CUTOFF]);
    tRise_m = Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]);
    tFall_m = Attributes::getReal(itsAttr[Attrib::Distribution::TFALL]);
    double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
    sigmaRise_m = tRise_m / tratio;
    sigmaFall_m = tFall_m / tratio;

    switch(distrTypeT_m) {
    case DistrTypeT::ASTRAFLATTOPTH: {
        double a = tPulseLengthFWHM_m / 2;
        double sig = tRise_m / 2;
        double inv_erf08 = 0.906193802436823; // erfinv(0.8)
        double sqr2 = sqrt(2.);
        double t = a - sqr2 * sig * inv_erf08;
        double tmps = sig;
        double tmpt = t;
        for (int i = 0; i < 10; ++ i) {
            sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
            t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
            sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
            tmps = sig;
            tmpt = t;
        }
        tEmission_m = tPulseLengthFWHM_m + 10 * sig;
        break;
    }
    case DistrTypeT::GUNGAUSSFLATTOPTH: {
        tEmission_m = tPulseLengthFWHM_m + (cutoff_m - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);
        break;
    }
    default:
        tEmission_m = 0.0;
    }
    return tEmission_m;
}

Inform &Distribution::printInfo(Inform &os) const {

    os << "\n"
       << "* ************* D I S T R I B U T I O N ********************************************" << endl;
    os << "* " << endl;
    if (OpalData::getInstance()->inRestartRun()) {
        os << "* In restart. Distribution read in from .h5 file." << endl;
    } else {
        if (addedDistributions_m.size() > 0)
            os << "* Main Distribution" << endl
               << "-----------------" << endl;

        if (particlesPerDist_m.empty())
            printDist(os, 0);
        else
            printDist(os, particlesPerDist_m.at(0));

        size_t distCount = 1;
        for (unsigned distIndex = 0; distIndex < addedDistributions_m.size(); distIndex++) {
            os << "* " << endl;
            os << "* Added Distribution #" << distCount << endl;
            os << "* ----------------------" << endl;
            addedDistributions_m.at(distIndex)->printDist(os, particlesPerDist_m.at(distCount));
            distCount++;
        }

        os << "* " << endl;
        if (numberOfEnergyBins_m > 0) {
            os << "* Number of energy bins    = " << numberOfEnergyBins_m << endl;

            //            if (numberOfEnergyBins_m > 1)
            //    printEnergyBins(os);
        }

        if (emitting_m) {
            os << "* Distribution is emitted. " << endl;
            os << "* Emission time            = " << tEmission_m << " [sec]" << endl;
            os << "* Time per bin             = " << tEmission_m / numberOfEnergyBins_m << " [sec]" << endl;
            os << "* Delta t during emission  = " << tBin_m / numberOfSampleBins_m << " [sec]" << endl;
            os << "* " << endl;
            printEmissionModel(os);
        } else {
            os << "* Distribution is injected." << endl;
        }
    }
    os << "* " << endl;
    os << "* **********************************************************************************" << endl;

    return os;
}

const PartData &Distribution::getReference() const {
    // Cast away const, to allow logically constant Distribution to update.
    const_cast<Distribution *>(this)->update();
    return particleRefData_m;
}

void Distribution::addDistributions() {
    /*
     * Move particle coordinates from added distributions to main distribution.
     */

    size_t idx = 1;
    std::vector<Distribution *>::iterator addedDistIt;
    for (addedDistIt = addedDistributions_m.begin();
         addedDistIt != addedDistributions_m.end(); ++ addedDistIt, ++ idx) {

        particlesPerDist_m[idx] = (*addedDistIt)->tOrZDist_m.size();

        for (double dist : (*addedDistIt)->getXDist()) {
            xDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseXDist();

        for (double dist : (*addedDistIt)->getBGxDist()) {
            pxDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseBGxDist();

        for (double dist : (*addedDistIt)->getYDist()) {
            yDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseYDist();

        for (double dist : (*addedDistIt)->getBGyDist()) {
            pyDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseBGyDist();

        for (double dist : (*addedDistIt)->getTOrZDist()) {
            tOrZDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseTOrZDist();

        for (double dist : (*addedDistIt)->getBGzDist()) {
            pzDist_m.push_back(dist);
        }
        (*addedDistIt)->eraseBGzDist();
    }
}

void Distribution::applyEmissionModel(double lowEnergyLimit, double &px, double &py, double &pz, std::vector<double> &additionalRNs) {

    switch (emissionModel_m) {

    case EmissionModelT::NONE:
        applyEmissModelNone(pz);
        break;
    case EmissionModelT::ASTRA:
        applyEmissModelAstra(px, py, pz, additionalRNs);
        break;
    case EmissionModelT::NONEQUIL:
        applyEmissModelNonEquil(lowEnergyLimit, px, py, pz, additionalRNs);
        break;
    default:
        break;
    }
}

void Distribution::applyEmissModelAstra(double &px, double &py, double &pz, std::vector<double> &additionalRNs) {

    double phi = 2.0 * acos(sqrt(additionalRNs[0]));
    double theta = 2.0 * Physics::pi * additionalRNs[1];

    px = pTotThermal_m * sin(phi) * cos(theta);
    py = pTotThermal_m * sin(phi) * sin(theta);
    pz = pTotThermal_m * std::abs(cos(phi));

}

void Distribution::applyEmissModelNone(double &pz) {
    pz += pTotThermal_m;
}

void Distribution::applyEmissModelNonEquil(double lowEnergyLimit,
                                           double &bgx,
                                           double &bgy,
                                           double &bgz,
                                           std::vector<double> &additionalRNs) {

    // Generate emission energy.
    double energy = 0.0;
    bool allow = false;

    const double expRelativeLaserEnergy = exp(laserEnergy_m / cathodeTemp_m);
    // double energyRange = emitEnergyUpperLimit_m - lowEnergyLimit;
    unsigned int counter = 0;
    while (!allow) {
        energy = lowEnergyLimit + additionalRNs[counter++] * emitEnergyUpperLimit_m;
        double randFuncValue = additionalRNs[counter++];
        double expRelativeEnergy = exp((energy - cathodeFermiEnergy_m) / cathodeTemp_m);
        double funcValue = ((1.0
                            - 1.0 / (1.0 + expRelativeEnergy * expRelativeLaserEnergy)) /
                            (1.0 + expRelativeEnergy));

        if (randFuncValue <= funcValue)
            allow = true;

        if (counter == additionalRNs.size()) {
            for (unsigned int i = 0; i < counter; ++ i) {
                additionalRNs[i] = gsl_rng_uniform(randGen_m);
            }

            counter = 0;
        }
    }

    while (additionalRNs.size() - counter < 2) {
        -- counter;
        additionalRNs[counter] = gsl_rng_uniform(randGen_m);
    }

    // Compute emission angles.
    double energyInternal = energy + laserEnergy_m;
    double energyExternal = energy - lowEnergyLimit; // uniformly distributed (?) value between 0 and emitEnergyUpperLimit_m

    double thetaInMax = acos(sqrt((lowEnergyLimit + laserEnergy_m) / (energyInternal)));
    double thetaIn = additionalRNs[counter++] * thetaInMax;
    double sinThetaOut = sin(thetaIn) * sqrt(energyInternal / energyExternal);
    double phi = Physics::two_pi * additionalRNs[counter];

    // Compute emission momenta.
    double betaGammaExternal
        = sqrt(pow(energyExternal / (Physics::m_e * 1.0e9) + 1.0, 2.0) - 1.0);

    bgx = betaGammaExternal * sinThetaOut * cos(phi);
    bgy = betaGammaExternal * sinThetaOut * sin(phi);
    bgz = betaGammaExternal * sqrt(1.0 - pow(sinThetaOut, 2.0));

}

void Distribution::calcPartPerDist(size_t numberOfParticles) {

    if (numberOfDistributions_m == 1) {
        particlesPerDist_m.push_back(numberOfParticles);
        return;
    }

    std::map<unsigned int, size_t> nPartFromFiles;
    double totalWeight = 0.0;
    for (unsigned int i = 0; i <= addedDistributions_m.size(); ++ i) {
        Distribution *currDist = this;
        if (i > 0)
            currDist = addedDistributions_m[i - 1];

        if (currDist->distrTypeT_m == DistrTypeT::FROMFILE) {
            std::ifstream inputFile;
            if (Ippl::myNode() == 0) {
                std::string fileName = Attributes::getString(currDist->itsAttr[Attrib::Distribution::FNAME]);
                inputFile.open(fileName.c_str());
            }
            size_t nPart = getNumberOfParticlesInFile(inputFile);
            nPartFromFiles.insert(std::make_pair(i, nPart));
            if (nPart > numberOfParticles) {
                throw OpalException("Distribution::calcPartPerDist",
                                    "Number of particles is too small");
            }

            numberOfParticles -= nPart;
        } else {
            totalWeight += currDist->getWeight();
        }
    }

    size_t numberOfCommittedPart = 0;
    for (unsigned int i = 0; i <= addedDistributions_m.size(); ++ i) {
        Distribution *currDist = this;
        if (i > 0)
            currDist = addedDistributions_m[i - 1];

        if (currDist->distrTypeT_m == DistrTypeT::FROMFILE) {
            particlesPerDist_m.push_back(nPartFromFiles[i]);
        } else {
            size_t particlesCurrentDist = numberOfParticles * currDist->getWeight() / totalWeight;
            particlesPerDist_m.push_back(particlesCurrentDist);
            numberOfCommittedPart += particlesCurrentDist;
        }
    }

    // Remaining particles go into first distribution that isn't FROMFILE
    if (numberOfParticles != numberOfCommittedPart) {
        size_t diffNumber = numberOfParticles - numberOfCommittedPart;
        for (unsigned int i = 0; i <= addedDistributions_m.size(); ++ i) {
            Distribution *currDist = this;
            if (i > 0)
                currDist = addedDistributions_m[i - 1];

            if (currDist->distrTypeT_m != DistrTypeT::FROMFILE) {
                particlesPerDist_m.at(i) += diffNumber;
                diffNumber = 0;
                break;
            }
        }
        if (diffNumber != 0) {
            throw OpalException("Distribution::calcPartPerDist",
                                "Particles can't be distributed to distributions in array");
        }
    }
}

void Distribution::checkEmissionParameters() {

    // There must be at least on energy bin for an emitted beam->
    numberOfEnergyBins_m
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[Attrib::Distribution::NBIN])));
    if (numberOfEnergyBins_m == 0)
        numberOfEnergyBins_m = 1;

    int emissionSteps = std::abs(static_cast<int> (Attributes::getReal(itsAttr[Attrib::Distribution::EMISSIONSTEPS])));

    // There must be at least one emission step.
    if (emissionSteps == 0)
        emissionSteps = 1;

    // Set number of sample bins per energy bin from the number of emission steps.
    numberOfSampleBins_m = static_cast<int> (std::ceil(emissionSteps / numberOfEnergyBins_m));
    if (numberOfSampleBins_m <= 0)
        numberOfSampleBins_m = 1;

    // Initialize emission counters.
    currentEnergyBin_m = 1;
    currentSampleBin_m = 0;

}

void Distribution::checkIfEmitted() {

    emitting_m = Attributes::getBool(itsAttr[Attrib::Distribution::EMITTED]);

    switch (distrTypeT_m) {

    case DistrTypeT::ASTRAFLATTOPTH:
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        emitting_m = true;
        break;
    default:
        break;
    }
}

void Distribution::checkParticleNumber(size_t &numberOfParticles) {

    size_t numberOfDistParticles = tOrZDist_m.size();
    reduce(numberOfDistParticles, numberOfDistParticles, OpAddAssign());

    if (numberOfDistParticles != numberOfParticles) {
        *gmsg << "\n--------------------------------------------------" << endl
              << "Warning!! The number of particles in the initial" << endl
              << "distribution is " << numberOfDistParticles << "." << endl << endl
              << "This is different from the number of particles" << endl
              << "defined by the BEAM command: " << numberOfParticles << endl << endl
              << "This is often happens when using a FROMFILE type" << endl
              << "distribution and not matching the number of" << endl
              << "particles in the particles file(s) with the number" << endl
              << "given in the BEAM command." << endl << endl
              << "The number of particles in the initial distribution" << endl
              << "(" << numberOfDistParticles << ") "
              << "will take precedence." << endl
              << "---------------------------------------------------\n" << endl;
    }
    numberOfParticles = numberOfDistParticles;
}

void Distribution::chooseInputMomentumUnits(InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits) {

    /*
     * Toggle what units to use for inputing momentum.
     */
    std::string inputUnits = Util::toUpper(Attributes::getString(itsAttr[Attrib::Distribution::INPUTMOUNITS]));
    if (inputUnits == "NONE")
        inputMoUnits_m = InputMomentumUnitsT::NONE;
    else if (inputUnits == "EV")
        inputMoUnits_m = InputMomentumUnitsT::EV;
    else
        inputMoUnits_m = inputMoUnits;

}

double Distribution::converteVToBetaGamma(double valueIneV, double massIneV) {
    double value = std::copysign(sqrt(pow(std::abs(valueIneV) / massIneV + 1.0, 2.0) - 1.0), valueIneV);
    if (std::abs(value) < std::numeric_limits<double>::epsilon())
        value = std::copysign(sqrt(2 * std::abs(valueIneV) / massIneV), valueIneV);

    return value;
}

void Distribution::createDistributionBinomial(size_t numberOfParticles, double massIneV) {

    setDistParametersBinomial(massIneV);
    generateBinomial(numberOfParticles);
}

void Distribution::createDistributionFlattop(size_t numberOfParticles, double massIneV) {

    setDistParametersFlattop(massIneV);

    if (emitting_m) {
        if (laserProfile_m == NULL)
            generateFlattopT(numberOfParticles);
        else
            generateFlattopLaserProfile(numberOfParticles);
    } else
        generateFlattopZ(numberOfParticles);
}

size_t Distribution::getNumberOfParticlesInFile(std::ifstream &inputFile) {

    size_t numberOfParticlesRead = 0;
    if (Ippl::myNode() == 0) {
        const boost::regex commentExpr("[[:space:]]*#.*");
        const std::string repl("");
        std::string line;
        std::stringstream linestream;
        long tempInt = 0;

        do {
            std::getline(inputFile, line);
            line = boost::regex_replace(line, commentExpr, repl);
        } while (line.length() == 0 && !inputFile.fail());

        linestream.str(line);
        linestream >> tempInt;
        if (tempInt <= 0) {
            throw OpalException("Distribution::getNumberOfParticlesInFile",
                                "The file '" +
                                Attributes::getString(itsAttr[Attrib::Distribution::FNAME]) +
                                "' does not seem to be an ASCII file containing a distribution.");
        }
        numberOfParticlesRead = static_cast<size_t>(tempInt);
    }
    reduce(numberOfParticlesRead, numberOfParticlesRead, OpAddAssign());

    return numberOfParticlesRead;
}

void Distribution::createDistributionFromFile(size_t numberOfParticles, double massIneV) {

    *gmsg << level3 << "\n"
          << "------------------------------------------------------------------------------------\n";
    *gmsg << "READ INITIAL DISTRIBUTION FROM FILE \""
          << Attributes::getString(itsAttr[Attrib::Distribution::FNAME])
          << "\"\n";
    *gmsg << "------------------------------------------------------------------------------------\n" << endl;

    // Data input file is only read by node 0.
    std::ifstream inputFile;
    std::string fileName = Attributes::getString(itsAttr[Attrib::Distribution::FNAME]);
    if (Ippl::myNode() == 0) {
        inputFile.open(fileName.c_str());
        if (inputFile.fail())
            throw OpalException("Distribution::createDistributionFromFile",
                                "Open file operation failed, please check if \""
                                + fileName
                                + "\" really exists.");
    }

    size_t numberOfParticlesRead = getNumberOfParticlesInFile(inputFile);
    /*
     * We read in the particle information using node zero, but save the particle
     * data to each processor in turn.
     */
    unsigned int saveProcessor = 0;

    pmean_m = 0.0;

    size_t numPartsToSend = 0;
    unsigned int distributeFrequency = 1000;
    size_t singleDataSize = (/*sizeof(int) +*/ 6 * sizeof(double));
    unsigned int dataSize = distributeFrequency * singleDataSize;
    std::vector<char> data;

    data.reserve(dataSize);

    const char* buffer;
    if (Ippl::myNode() == 0) {
        char lineBuffer[1024];
        unsigned int numParts = 0;
        while (!inputFile.eof()) {
            inputFile.getline(lineBuffer, 1024);

            Vector_t R(0.0), P(0.0);

            std::istringstream line(lineBuffer);
            line >> R(0);
            if (line.rdstate()) break;
            line >> P(0);
            line >> R(1);
            line >> P(1);
            line >> R(2);
            line >> P(2);

            if (saveProcessor >= (unsigned)Ippl::getNodes())
                saveProcessor = 0;

            if (inputMoUnits_m == InputMomentumUnitsT::EV) {
                P(0) = converteVToBetaGamma(P(0), massIneV);
                P(1) = converteVToBetaGamma(P(1), massIneV);
                P(2) = converteVToBetaGamma(P(2), massIneV);
            }

            pmean_m += P;

            if (saveProcessor > 0u) {
                buffer = reinterpret_cast<const char*>(&R[0]);
                data.insert(data.end(), buffer, buffer + 3 * sizeof(double));
                buffer = reinterpret_cast<const char*>(&P[0]);
                data.insert(data.end(), buffer, buffer + 3 * sizeof(double));
                ++ numPartsToSend;

                if (numPartsToSend % distributeFrequency == 0) {
                    MPI_Bcast(&dataSize, 1, MPI_INT, 0, Ippl::getComm());
                    MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());
                    numPartsToSend = 0;

                    std::vector<char>().swap(data);
                    data.reserve(dataSize);
                }
            } else {
                xDist_m.push_back(R(0));
                yDist_m.push_back(R(1));
                tOrZDist_m.push_back(R(2));
                pxDist_m.push_back(P(0));
                pyDist_m.push_back(P(1));
                pzDist_m.push_back(P(2));
            }

            ++ numParts;
            ++ saveProcessor;
        }

        dataSize = (numberOfParticlesRead == numParts? data.size(): std::numeric_limits<unsigned int>::max());
        MPI_Bcast(&dataSize, 1, MPI_INT, 0, Ippl::getComm());
        if (numberOfParticlesRead != numParts) {
            throw OpalException("Distribution::createDistributionFromFile",
                                "Found " +
                                std::to_string(numParts) +
                                " particles in file '" +
                                fileName +
                                "' instead of " +
                                std::to_string(numberOfParticlesRead));
        }
        MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());

    } else {
        do {
            MPI_Bcast(&dataSize, 1, MPI_INT, 0, Ippl::getComm());
            if (dataSize == std::numeric_limits<unsigned int>::max()) {
                throw OpalException("Distribution::createDistributionFromFile",
                                    "Couldn't find " +
                                    std::to_string(numberOfParticlesRead) +
                                    " particles in file '" +
                                    fileName + "'");
            }
            MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());

            size_t i = 0;
            while (i < dataSize) {

                if (saveProcessor + 1 == (unsigned) Ippl::myNode()) {
                    const double *tmp = reinterpret_cast<const double*>(&data[i]);
                    xDist_m.push_back(tmp[0]);
                    yDist_m.push_back(tmp[1]);
                    tOrZDist_m.push_back(tmp[2]);
                    pxDist_m.push_back(tmp[3]);
                    pyDist_m.push_back(tmp[4]);
                    pzDist_m.push_back(tmp[5]);
                    i += 6 * sizeof(double);
                } else {
                    i += singleDataSize;
                }

                ++ saveProcessor;
                if (saveProcessor + 1 >= (unsigned) Ippl::getNodes()) {
                    saveProcessor = 0;
                }
            }
        } while (dataSize == distributeFrequency * singleDataSize);
    }

    pmean_m /= numberOfParticlesRead;
    reduce(pmean_m, pmean_m, OpAddAssign());

    if (Ippl::myNode() == 0)
        inputFile.close();
}


void Distribution::createMatchedGaussDistribution(size_t numberOfParticles, double massIneV) {

    /*
      ToDo:
      - eliminate physics and error
    */

    std::string LineName = Attributes::getString(itsAttr[Attrib::Distribution::LINE]);
    if (LineName == "") return;

    const BeamSequence* LineSequence = BeamSequence::find(LineName);

    if (LineSequence == NULL)
        throw OpalException("Distribution::CreateMatchedGaussDistribution",
                            "didn't find any Cyclotron element in line");

    SpecificElementVisitor<Cyclotron> CyclotronVisitor(*LineSequence->fetchLine());
    CyclotronVisitor.execute();
    size_t NumberOfCyclotrons = CyclotronVisitor.size();

    if (NumberOfCyclotrons > 1) {
        throw OpalException("Distribution::createMatchedGaussDistribution",
                            "I am confused, found more than one Cyclotron element in line");
    }
    if (NumberOfCyclotrons == 0) {
        throw OpalException("Distribution::createMatchedGaussDistribution",
                            "didn't find any Cyclotron element in line");
    }

    /* FIXME we need to remove const-ness otherwise we can't update the injection radius
     * and injection radial momentum. However, there should be a better solution ..
     */
    Cyclotron* CyclotronElement = const_cast<Cyclotron*>(CyclotronVisitor.front());

    bool full    = !Attributes::getBool(itsAttr[Attrib::Distribution::SECTOR]);
    int Nint     = (int)Attributes::getReal(itsAttr[Attrib::Distribution::NSTEPS]);
    int Nsectors = (int)Attributes::getReal(itsAttr[Attrib::Distribution::NSECTORS]);

    if ( Nint < 0 )
        throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                            "Negative number of integration steps");

    if ( Nsectors < 0 )
        throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                            "Negative number of sectors");

    if ( Nsectors > 1 && full == false )
        throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                            "Averaging over sectors can only be done with SECTOR=FALSE");

    *gmsg << "* ----------------------------------------------------" << endl;
    *gmsg << "* About to find closed orbit and matched distribution " << endl;
    *gmsg << "* I= " << I_m*1E3 << " (mA)  E= " << E_m*1E-6 << " (MeV)" << endl;
    *gmsg << "* EX= " << Attributes::getReal(itsAttr[Attrib::Distribution::EX])
          << "  EY= " << Attributes::getReal(itsAttr[Attrib::Distribution::EY])
          << "  ET= " << Attributes::getReal(itsAttr[Attrib::Distribution::ET]) << endl;
    if ( full ) {
        *gmsg << "* SECTOR: " << "match using all sectors, with" << endl
              << "* NSECTORS = " << Nsectors << " to average the field over" << endl;
    }
    else
        *gmsg << "* SECTOR: " << "match using single sector" << endl;

    *gmsg << "* NSTEPS = " << Nint << endl
          << "* HN= "      << CyclotronElement->getCyclHarm()
          << "  PHIINIT= " << CyclotronElement->getPHIinit()  << endl
          << "* ----------------------------------------------------" << endl;

    if ( CyclotronElement->getFMLowE()  < 0 ||
         CyclotronElement->getFMHighE() < 0 )
    {
        throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                            "Missing attributes 'FMLOWE' and/or 'FMHIGHE' in "
                            "'CYCLOTRON' definition.");
    }

    std::size_t maxitCOF =
        Attributes::getReal(itsAttr[Attrib::Distribution::MAXSTEPSCO]);

    double rguess =
        Attributes::getReal(itsAttr[Attrib::Distribution::RGUESS]);
    
    double denergy = 1000.0 *
        Attributes::getReal(itsAttr[Attrib::Distribution::DENERGY]);
    
    if ( denergy < 0.0 )
        throw OpalException("Distribution:CreateMatchedGaussDistribution()",
                            "DENERGY < 0");

    double accuracy =
        Attributes::getReal(itsAttr[Attrib::Distribution::RESIDUUM]);

    if ( Options::cloTuneOnly ) {
        *gmsg << "* Stopping after closed orbit and tune calculation" << endl;
        typedef std::vector<double> container_t;
        typedef boost::numeric::odeint::runge_kutta4<container_t> rk4_t;
        typedef ClosedOrbitFinder<double,unsigned int, rk4_t> cof_t;

        cof_t cof(massIneV*1E-6, Nint, CyclotronElement, false, Nsectors);
        cof.findOrbit(accuracy, maxitCOF, E_m*1E-6, denergy, rguess, true);

        std::exit(0);
    }

    bool writeMap = true;

    typedef SigmaGenerator<double, unsigned int> sGenerator_t;
    sGenerator_t *siggen = new sGenerator_t(I_m,
                                            Attributes::getReal(itsAttr[Attrib::Distribution::EX])*1E6,
                                            Attributes::getReal(itsAttr[Attrib::Distribution::EY])*1E6,
                                            Attributes::getReal(itsAttr[Attrib::Distribution::ET])*1E6,
                                            E_m*1E-6,
                                            massIneV*1E-6,
                                            CyclotronElement,
                                            Nint,
                                            Attributes::getReal(itsAttr[Attrib::Distribution::ORDERMAPS]),
                                            writeMap);

    if (siggen->match(accuracy,
                      Attributes::getReal(itsAttr[Attrib::Distribution::MAXSTEPSSI]),
                      maxitCOF,
                      CyclotronElement,
                      denergy,
                      rguess,
                      false, full))  {

        std::array<double,3> Emit = siggen->getEmittances();

        if (rguess > 0)
            *gmsg << "* RGUESS " << rguess/1000.0 << " (m) " << endl;

        *gmsg << "* Converged (Ex, Ey, Ez) = (" << Emit[0] << ", " << Emit[1] << ", "
              << Emit[2] << ") pi mm mrad for E= " << E_m*1E-6 << " (MeV)" << endl;
        *gmsg << "* Sigma-Matrix " << endl;

        for (unsigned int i = 0; i < siggen->getSigma().size1(); ++ i) {
            *gmsg << std::setprecision(4)  << std::setw(11) << siggen->getSigma()(i,0);
            for (unsigned int j = 1; j < siggen->getSigma().size2(); ++ j) {
                if (siggen->getSigma()(i,j) < 10e-12){
                    *gmsg << " & " <<  std::setprecision(4)  << std::setw(11) << 0.0;
                }
                else{
                   *gmsg << " & " <<  std::setprecision(4)  << std::setw(11) << siggen->getSigma()(i,j);
                }

            }
            *gmsg << " \\\\" << endl;
        }

        /*

          Now setup the distribution generator
          Units of the Sigma Matrix:

          spatial: mm
          moment:  rad

        */
        double gamma = E_m / massIneV + 1.0;
        double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));

        auto sigma = siggen->getSigma();
        // change units from mm to m
        for (unsigned int i = 0; i < 6; ++ i)
            for (unsigned int j = 0; j < 6; ++ j) sigma(i, j) *= 1e-6;

        for (unsigned int i = 0; i < 3; ++ i) {
            if ( sigma(2 * i, 2 * i) < 0 || sigma(2 * i + 1, 2 * i + 1) < 0 )
                throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                                    "Negative value on the diagonal of the sigma matrix.");
        }

        sigmaR_m[0] = std::sqrt(sigma(0, 0));
        sigmaP_m[0] = std::sqrt(sigma(1, 1))*beta*gamma;
        sigmaR_m[2] = std::sqrt(sigma(2, 2));
        sigmaP_m[2] = std::sqrt(sigma(3, 3))*beta*gamma;
        sigmaR_m[1] = std::sqrt(sigma(4, 4));

        //p_l^2 = [(delta+1)*beta*gamma]^2 - px^2 - pz^2

        double pl2 = (std::sqrt(sigma(5,5)) + 1)*(std::sqrt(sigma(5,5)) + 1)*beta*gamma*beta*gamma -
                      sigmaP_m[0]*sigmaP_m[0] - sigmaP_m[2]*sigmaP_m[2];

        double pl = std::sqrt(pl2);
        sigmaP_m[1] = gamma*(pl - beta*gamma);

        correlationMatrix_m(1, 0) = sigma(0, 1) / (sqrt(sigma(0, 0) * sigma(1, 1)));
        correlationMatrix_m(3, 2) = sigma(2, 3) / (sqrt(sigma(2, 2) * sigma(3, 3)));
        correlationMatrix_m(5, 4) = sigma(4, 5) / (sqrt(sigma(4, 4) * sigma(5, 5)));
        correlationMatrix_m(4, 0) = sigma(0, 4) / (sqrt(sigma(0, 0) * sigma(4, 4)));
        correlationMatrix_m(4, 1) = sigma(1, 4) / (sqrt(sigma(1, 1) * sigma(4, 4)));
        correlationMatrix_m(5, 0) = sigma(0, 5) / (sqrt(sigma(0, 0) * sigma(5, 5)));
        correlationMatrix_m(5, 1) = sigma(1, 5) / (sqrt(sigma(1, 1) * sigma(5, 5)));

        createDistributionGauss(numberOfParticles, massIneV);

        // update injection radius and radial momentum
        CyclotronElement->setRinit(siggen->getInjectionRadius() * 1.0e3);
        CyclotronElement->setPRinit(siggen->getInjectionMomentum());
    }
    else {
        *gmsg << "* Not converged for " << E_m*1E-6 << " MeV" << endl;

        delete siggen;

        throw OpalException("Distribution::CreateMatchedGaussDistribution",
                            "didn't find any matched distribution.");
    }

    delete siggen;
}

void Distribution::createDistributionGauss(size_t numberOfParticles, double massIneV) {

    setDistParametersGauss(massIneV);

    if (emitting_m) {
        generateTransverseGauss(numberOfParticles);
        generateLongFlattopT(numberOfParticles);
    } else {
        generateGaussZ(numberOfParticles);
    }
}

void  Distribution::createBoundaryGeometry(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg) {

    int pc = 0;
    size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
    size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
    if (Ippl::myNode() == 0)
        N_mean += N_extra;
    size_t count = 0;
    size_t lowMark = beam->getLocalNum();
    if (bg.getN() != 0) {

        for (size_t i = 0; i < bg.getN(); i++) {
            if (pc == Ippl::myNode()) {
                if (count < N_mean) {
                    beam->create(1);
                    if (pc != 0)
                        beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                    else
                        beam->R[lowMark + count] = bg.getCooridinate(count);
                    beam->P[lowMark + count] = Vector_t(0.0);
                    beam->Bin[lowMark + count] = 0;
                    beam->PType[lowMark + count] = ParticleType::FIELDEMISSION;
                    beam->TriID[lowMark + count] = 0;
                    beam->Q[lowMark + count] = beam->getChargePerParticle();
                    beam->Ef[lowMark + count] = Vector_t(0.0);
                    beam->Bf[lowMark + count] = Vector_t(0.0);
                    beam->dt[lowMark + count] = beam->getdT();
                    count ++;
                }
            }
            pc++;
            if (pc == Ippl::getNodes())
                pc = 0;
        }
    }
    bg.clearCooridinateArray();
    beam->boundp();
    *gmsg << &beam << endl;
}

void Distribution::createOpalCycl(PartBunchBase<double, 3> *beam,
                                  size_t numberOfParticles,
                                  double current, const Beamline &bl) {

    /*
     *  setup data for matched distribution generation
     */

    E_m = (beam->getInitialGamma()-1.0)*beam->getM();
    I_m = current;

    /*
      Fixme:

      avrgpz_m = beam->getP()/beam->getM();
    */
    size_t numberOfPartToCreate = numberOfParticles;
    totalNumberParticles_m = numberOfParticles;
    if (beam->getTotalNum() != 0) {
        numberOfPartToCreate = beam->getLocalNum();
    }

    // Setup particle bin structure.
    setupParticleBins(beam->getM(),beam);

    /*
     * Set what units to use for input momentum units. Default is
     * eV.
     */
    chooseInputMomentumUnits(InputMomentumUnitsT::EV);

    /*
     * Determine the number of particles for each distribution. For OPAL-cycl
     * there are currently no arrays of distributions supported
     */
    calcPartPerDist(numberOfParticles);

    // Set distribution type.
    setDistType();

    // Emitting particles in not supported in OpalCyclT.
    emitting_m = false;

    // Create distribution.
    create(numberOfPartToCreate, beam->getM());

    // this variable is needed to be compatible with OPAL-T
    particlesPerDist_m.push_back(tOrZDist_m.size());

    shiftDistCoordinates(beam->getM());

    // Setup energy bins.
    if (numberOfEnergyBins_m > 0) {
        fillParticleBins();
        beam->setPBins(energyBins_m);
    }

    // Check number of particles in distribution.
    checkParticleNumber(numberOfParticles);

    initializeBeam(beam);
    writeOutFileHeader();

    injectBeam(beam);

    OpalData::getInstance()->addProblemCharacteristicValue("NP", numberOfParticles);
}

void Distribution::createOpalE(Beam *beam,
                               std::vector<Distribution *> addedDistributions,
                               EnvelopeBunch *envelopeBunch,
                               double distCenter,
                               double Bz0) {

    IpplTimings::startTimer(envelopeBunch->distrCreate_m);

    double beamEnergy = beam->getMass() * (beam->getGamma() - 1.0) * 1.0e9;
    numberOfEnergyBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[Attrib::Distribution::NBIN])));

    /*
     * Set what units to use for input momentum units. Default is
     * unitless (i.e. BetaXGamma, BetaYGamma, BetaZGamma).
     */
    chooseInputMomentumUnits(InputMomentumUnitsT::NONE);

    setDistType();

    // Check if this is to be an emitted beam->
    checkIfEmitted();

    switch (distrTypeT_m) {

    case DistrTypeT::FLATTOP:
        setDistParametersFlattop(beam->getMass());
        beamEnergy = Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]);
        break;
    case DistrTypeT::GAUSS:
        setDistParametersGauss(beam->getMass());
        break;
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        setDistParametersFlattop(beam->getMass());
        beamEnergy = Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]);
        break;
    default:
        *gmsg << "Only FLATTOP, GAUSS and GUNGAUSSFLATTOPTH distribution types supported " << endl
              << "in envelope mode. Assuming FLATTOP." << endl;
        distrTypeT_m = DistrTypeT::FLATTOP;
        setDistParametersFlattop(beam->getMass());
        break;
    }

    tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
        * (sigmaTRise_m + sigmaTFall_m);
    double beamWidth = tEmission_m * Physics::c * sqrt(1.0 - (1.0 / pow(beam->getGamma(), 2)));
    double beamCenter = -1.0 * beamWidth / 2.0;

    envelopeBunch->initialize(beam->getNumberOfSlices(),
                              beam->getCharge(),
                              beamEnergy,
                              beamWidth,
                              tEmission_m,
                              0.9,
                              beam->getCurrent(),
                              beamCenter,
                              sigmaR_m[0],
                              sigmaR_m[1],
                              0.0,
                              0.0,
                              Bz0,
                              numberOfEnergyBins_m);

    IpplTimings::stopTimer(envelopeBunch->distrCreate_m);
}

void Distribution::createOpalT(PartBunchBase<double, 3> *beam,
                               std::vector<Distribution *> addedDistributions,
                               size_t &numberOfParticles) {

    addedDistributions_m = addedDistributions;
    createOpalT(beam, numberOfParticles);
}

void Distribution::createOpalT(PartBunchBase<double, 3> *beam,
                               size_t &numberOfParticles) {

    IpplTimings::startTimer(beam->distrCreate_m);

    // This is PC from BEAM
    double deltaP = Attributes::getReal(itsAttr[Attrib::Distribution::OFFSETP]);
    if (inputMoUnits_m == InputMomentumUnitsT::EV) {
        deltaP = converteVToBetaGamma(deltaP, beam->getM());
    }

    avrgpz_m = beam->getP()/beam->getM() + deltaP;

    totalNumberParticles_m = numberOfParticles;

    /*
     * Set what units to use for input momentum units. Default is
     * unitless (i.e. BetaXGamma, BetaYGamma, BetaZGamma).
     */
    chooseInputMomentumUnits(InputMomentumUnitsT::NONE);

    // Set distribution type(s).
    setDistType();
    for (Distribution* addedDist : addedDistributions_m)
        addedDist->setDistType();

    /*
     * Determine the number of particles for each distribution. Note
     * that if a distribution is generated from an input file, then
     * the number of particles in that file will override what is
     * determined here.
     */
    calcPartPerDist(numberOfParticles);

    // Check if this is to be an emitted beam->
    checkIfEmitted();

    /*
     * Force added distributions to the same emission condition as the main
     * distribution.
     */
    for (Distribution* addedDist : addedDistributions_m)
        addedDist->setDistToEmitted(emitting_m);

    if (emitting_m)
        setupEmissionModel(beam);

    // Create distributions.
    create(particlesPerDist_m.at(0), beam->getM());

    size_t distCount = 1;
    for (Distribution* addedDist : addedDistributions_m) {
        addedDist->create(particlesPerDist_m.at(distCount), beam->getM());
        distCount++;
    }

    // Move added distribution particles to main distribution.
    addDistributions();

    if (emitting_m && emissionModel_m == EmissionModelT::NONE)
        setupEmissionModelNone(beam);

    // Check number of particles in distribution.
    checkParticleNumber(numberOfParticles);

    if (emitting_m) {
        checkEmissionParameters();
    } else {
        if (distrTypeT_m != DistrTypeT::FROMFILE) {
            pmean_m = Vector_t(0, 0, avrgpz_m);
        }
    }

    /*
     * Find max. and min. particle positions in the bunch. For the
     * case of an emitted beam these will be in time (seconds). For
     * an injected beam in z (meters).
     */
    double maxTOrZ = getMaxTOrZ();
    double minTOrZ = getMinTOrZ();

    /*
     * Set emission time and/or shift distributions if necessary.
     * For an emitted beam we place all particles at negative time.
     * For an injected beam we just ensure that there are no
     * particles at z < 0.
     */

    if (emitting_m) {
        setEmissionTime(maxTOrZ, minTOrZ);
    }
    shiftBeam(maxTOrZ, minTOrZ);

    shiftDistCoordinates(beam->getM());

    if (numberOfEnergyBins_m > 0) {
        setupEnergyBins(maxTOrZ, minTOrZ);
        fillEBinHistogram();
    }

    initializeBeam(beam);
    writeOutFileHeader();

    if (emitting_m && Options::cZero) {
        std::vector<std::vector<double> > mirrored;
        const auto end = additionalRNs_m.end();

        if (emissionModel_m == EmissionModelT::ASTRA ||
            distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH ||
            distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH) {

            for (auto it = additionalRNs_m.begin(); it != end; ++ it) {
                std::vector<double> tmp;
                tmp.push_back((*it).front());
                tmp.push_back((*it).back() + 0.5);
                mirrored.push_back(tmp);
            }
        } else {
            size_t numAdditionals = additionalRNs_m.front().size() / 2;
            for (auto it = additionalRNs_m.begin(); it != end; ++ it) {
                std::vector<double> tmp((*it).begin() + numAdditionals, (*it).end());
                mirrored.push_back(tmp);
                (*it).erase((*it).begin() + numAdditionals, (*it).end());
            }
        }

        additionalRNs_m.insert(additionalRNs_m.end(), mirrored.begin(), mirrored.end());
    }

    /*
     * If this is an injected beam, we create particles right away.
     * Emitted beams get created during the course of the simulation.
     */
    if (!emitting_m)
        injectBeam(beam);

    OpalData::getInstance()->addProblemCharacteristicValue("NP", numberOfParticles);
    IpplTimings::stopTimer(beam->distrCreate_m);
}

/**
 * Here we emit particles from the cathode.

 A typical integration time step, \f$\Delta t\f$, is broken down into 3 sub-steps:

 1) Drift particles for \f$\frac{\Delta t}{2}\f$.

 2) Calculate fields and advance momentum.

 3) Drift particles for \f$\frac{\Delta t}{2}\f$ at the new momentum to complete the
 full time step.

 The difficulty for emission is that at the cathode position there is a step function discontinuity in the  fields. If we
 apply the typical integration time step across this boundary, we get an artificial numerical bunching of the beam, especially
 at very high accelerating fields. This function takes the cathode position boundary into account in order to achieve
 smoother particle emission.

 During an emission step, an integral number of time bins from the distribution histogram are emitted. However, each particle
 contained in those time bins will actually be emitted from the cathode at a different time, so will only spend some fraction
 of the total time step, \f$\Delta t_{full-timestep}\f$, in the simulation. The trick to emission is to give each particle
 a unique time step, \f$Delta t_{temp}\f$, that is equal to the actual time during the emission step that the particle
 exists in the simulation. For the next integration time step, the particle's time step is set back to the global time step,
 \f$\Delta t_{full-timestep}\f$.
*/
size_t Distribution::emitParticles(PartBunchBase<double, 3> *beam, double eZ) {

    // Number of particles that have already been emitted and are on this processor.
    size_t numberOfEmittedParticles = beam->getLocalNum();
    size_t oldNumberOfEmittedParticles = numberOfEmittedParticles;

    if (!tOrZDist_m.empty() && emitting_m) {

        /*
         * Loop through emission beam coordinate containers and emit particles with
         * the appropriate time coordinate. Once emitted, remove particle from the
         * "to be emitted" list.
         */
        std::vector<size_t> particlesToBeErased;
        double phiEffective = (cathodeWorkFunc_m
                                     - sqrt(std::max(0.0, (Physics::q_e * beam->getQ() * eZ) /
                                                     (4.0 * Physics::pi * Physics::epsilon_0))));
        double lowEnergyLimit = cathodeFermiEnergy_m + phiEffective - laserEnergy_m;

        for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {

            // Advance particle time.
            tOrZDist_m.at(particleIndex) += beam->getdT();

            // If particle time is now greater than zero, we emit it.
            if (tOrZDist_m.at(particleIndex) >= 0.0) {

                particlesToBeErased.push_back(particleIndex);

                beam->create(1);
                double deltaT = tOrZDist_m.at(particleIndex);
                beam->dt[numberOfEmittedParticles] = deltaT;

                double oneOverCDt = 1.0 / (Physics::c * deltaT);

                double px = pxDist_m.at(particleIndex);
                double py = pyDist_m.at(particleIndex);
                double pz = pzDist_m.at(particleIndex);
                std::vector<double> additionalRNs;
                if (additionalRNs_m.size() > particleIndex) {
                    additionalRNs = additionalRNs_m[particleIndex];
                } else {
                    throw OpalException("Distribution::emitParticles",
                                        "not enough additional particles");
                }
                applyEmissionModel(lowEnergyLimit, px, py, pz, additionalRNs);

                double particleGamma
                    = std::sqrt(1.0
                                + std::pow(px, 2.0)
                                + std::pow(py, 2.0)
                                + std::pow(pz, 2.0));

                beam->R[numberOfEmittedParticles]
                    = Vector_t(oneOverCDt * (xDist_m.at(particleIndex)
                                             + px * deltaT * Physics::c / (2.0 * particleGamma)),
                               oneOverCDt * (yDist_m.at(particleIndex)
                                             + py * deltaT * Physics::c / (2.0 * particleGamma)),
                               pz / (2.0 * particleGamma));
                beam->P[numberOfEmittedParticles]
                    = Vector_t(px, py, pz);
                beam->Bin[numberOfEmittedParticles] = currentEnergyBin_m - 1;
                beam->Q[numberOfEmittedParticles] = beam->getChargePerParticle();
                beam->Ef[numberOfEmittedParticles] = Vector_t(0.0);
                beam->Bf[numberOfEmittedParticles] = Vector_t(0.0);
                beam->PType[numberOfEmittedParticles] = ParticleType::REGULAR;
                beam->TriID[numberOfEmittedParticles] = 0;
                numberOfEmittedParticles++;

                beam->iterateEmittedBin(currentEnergyBin_m - 1);

                // Save particles to vectors for writing initial distribution.
                xWrite_m.push_back(xDist_m.at(particleIndex));
                pxWrite_m.push_back(px);
                yWrite_m.push_back(yDist_m.at(particleIndex));
                pyWrite_m.push_back(py);
                tOrZWrite_m.push_back(-(beam->getdT() - deltaT + currentEmissionTime_m));
                pzWrite_m.push_back(pz);
                binWrite_m.push_back(currentEnergyBin_m);
            }
        }

        // Now erase particles that were emitted.
        std::vector<size_t>::reverse_iterator ptbErasedIt;
        for (ptbErasedIt = particlesToBeErased.rbegin();
             ptbErasedIt < particlesToBeErased.rend();
             ++ptbErasedIt) {

            /*
             * We don't use the vector class function erase because it
             * is much slower than doing a swap and popping off the end
             * of the vector.
             */
            std::swap(xDist_m.at(*ptbErasedIt), xDist_m.back());
            std::swap(pxDist_m.at(*ptbErasedIt), pxDist_m.back());
            std::swap(yDist_m.at(*ptbErasedIt), yDist_m.back());
            std::swap(pyDist_m.at(*ptbErasedIt), pyDist_m.back());
            std::swap(tOrZDist_m.at(*ptbErasedIt), tOrZDist_m.back());
            std::swap(pzDist_m.at(*ptbErasedIt), pzDist_m.back());
            if (additionalRNs_m.size() == xDist_m.size()) {
                std::swap(additionalRNs_m.at(*ptbErasedIt), additionalRNs_m.back());

                additionalRNs_m.pop_back();
            }

            xDist_m.pop_back();
            pxDist_m.pop_back();
            yDist_m.pop_back();
            pyDist_m.pop_back();
            tOrZDist_m.pop_back();
            pzDist_m.pop_back();

        }

        /*
         * Set energy bin to emitted if all particles in the bin have been emitted.
         * However, be careful with the last energy bin. We cannot emit it until all
         * of the particles have been accounted for. So when on the last bin, keep it
         * open for the rest of the beam->
         */
        currentSampleBin_m++;
        if (currentSampleBin_m == numberOfSampleBins_m) {

            INFOMSG(level3 << "* Bin number: "
                    << currentEnergyBin_m
                    << " has emitted all particles (new emit)." << endl);
            currentSampleBin_m = 0;
            currentEnergyBin_m++;

        }

        /*
         * Set beam to emitted. Make sure temporary storage is cleared.
         */
        if (currentEnergyBin_m > numberOfEnergyBins_m || tOrZDist_m.empty()) {
            emitting_m = false;

            xDist_m.clear();
            pxDist_m.clear();
            yDist_m.clear();
            pyDist_m.clear();
            tOrZDist_m.clear();
            pzDist_m.clear();

            currentEnergyBin_m = numberOfEnergyBins_m;
        }

    }
    currentEmissionTime_m += beam->getdT();

    // Write emitted particles to file.
    writeOutFileEmission();

    size_t currentEmittedParticles = numberOfEmittedParticles - oldNumberOfEmittedParticles;
    reduce(currentEmittedParticles, currentEmittedParticles, OpAddAssign());
    totalNumberEmittedParticles_m += currentEmittedParticles;

    return currentEmittedParticles;

}

void Distribution::eraseXDist() {
    xDist_m.erase(xDist_m.begin(), xDist_m.end());
}

void Distribution::eraseBGxDist() {
    pxDist_m.erase(pxDist_m.begin(), pxDist_m.end());
}

void Distribution::eraseYDist() {
    yDist_m.erase(yDist_m.begin(), yDist_m.end());
}

void Distribution::eraseBGyDist() {
    pyDist_m.erase(pyDist_m.begin(), pyDist_m.end());
}

void Distribution::eraseTOrZDist() {
    tOrZDist_m.erase(tOrZDist_m.begin(), tOrZDist_m.end());
}

void Distribution::eraseBGzDist() {
    pzDist_m.erase(pzDist_m.begin(), pzDist_m.end());
}

void Distribution::fillEBinHistogram() {
    double upper = 0.0;
    double lower = 0.0;
    gsl_histogram_get_range(energyBinHist_m,
                            gsl_histogram_bins(energyBinHist_m) - 1,
                            &lower, &upper);
    const size_t numberOfParticles = tOrZDist_m.size();
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double &tOrZ = tOrZDist_m.at(partIndex);

        if (gsl_histogram_increment(energyBinHist_m, tOrZ) == GSL_EDOM) {
            gsl_histogram_increment(energyBinHist_m, 0.5 * (lower + upper));
        }
    }

    reduce(energyBinHist_m->bin, energyBinHist_m->bin + energyBinHist_m->n, energyBinHist_m->bin, OpAddAssign());
}

void Distribution::fillParticleBins() {

    for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {

        std::vector<double> partCoord;

        partCoord.push_back(xDist_m.at(particleIndex));
        partCoord.push_back(yDist_m.at(particleIndex));
        partCoord.push_back(tOrZDist_m.at(particleIndex));
        partCoord.push_back(pxDist_m.at(particleIndex));
        partCoord.push_back(pyDist_m.at(particleIndex));
        partCoord.push_back(pzDist_m.at(particleIndex));
        partCoord.push_back(0.0);

        energyBins_m->fill(partCoord);

    }

    energyBins_m->sortArray();
}

size_t Distribution::findEBin(double tOrZ) {

    if (tOrZ <= gsl_histogram_min(energyBinHist_m)) {
        return 0;
    } else if (tOrZ >= gsl_histogram_max(energyBinHist_m)) {
        return numberOfEnergyBins_m - 1;
    } else {
        size_t binNumber;
        gsl_histogram_find(energyBinHist_m, tOrZ, &binNumber);
        return binNumber / numberOfSampleBins_m;
    }
}

void Distribution::generateAstraFlattopT(size_t numberOfParticles) {

    /*
     * Legacy function to support the ASTRAFLATTOPOTH
     * distribution type.
     */
    checkEmissionParameters();

    gsl_qrng *quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, 2);

    int numberOfSampleBins
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[Attrib::Legacy::Distribution::SBIN])));
    int numberOfEnergyBins
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[Attrib::Distribution::NBIN])));

    int binTotal = numberOfSampleBins * numberOfEnergyBins;

    double *distributionTable = new double[binTotal + 1];

    double a = tPulseLengthFWHM_m / 2.;
    double sig = tRise_m / 2.;
    double inv_erf08 = 0.906193802436823; // erfinv(0.8)
    double sqr2 = sqrt(2.0);
    double t = a - sqr2 * sig * inv_erf08;
    double tmps = sig;
    double tmpt = t;

    for (int i = 0; i < 10; ++ i) {
        sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
        t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
        sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
        tmps = sig;
        tmpt = t;
    }

    /*
     * Emission time is set here during distribution particle creation only for
     * the ASTRAFLATTOPTH distribution type.
     */
    tEmission_m = tPulseLengthFWHM_m + 10. * sig;
    tBin_m = tEmission_m / numberOfEnergyBins;

    double lo = -tBin_m / 2.0 * numberOfEnergyBins;
    double hi = tBin_m / 2.0 * numberOfEnergyBins;
    double dx = tBin_m / numberOfSampleBins;
    double x = lo;
    double tot = 0;
    double weight = 2.0;

    // sample the function that describes the histogram of the requested distribution
    for (int i = 0; i < binTotal + 1; ++ i, x += dx, weight = 6. - weight) {
        distributionTable[i] = gsl_sf_erf((x + a) / (sqrt(2) * sig)) - gsl_sf_erf((x - a) / (sqrt(2) * sig));
        tot += distributionTable[i] * weight;
    }
    tot -= distributionTable[binTotal] * (5. - weight);
    tot -= distributionTable[0];

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);
    double tCoord = 0.0;

    int effectiveNumParticles = 0;
    int largestBin = 0;
    std::vector<int> numParticlesInBin(numberOfEnergyBins,0);
    for (int k = 0; k < numberOfEnergyBins; ++ k) {
        double loc_fraction = -distributionTable[numberOfSampleBins * k] / tot;

        weight = 2.0;
        for (int i = numberOfSampleBins * k; i < numberOfSampleBins * (k + 1) + 1;
            ++ i, weight = 6. - weight) {
            loc_fraction += distributionTable[i] * weight / tot;
        }
        loc_fraction -= distributionTable[numberOfSampleBins * (k + 1)]
            * (5. - weight) / tot;
        numParticlesInBin[k] = static_cast<int>(std::floor(loc_fraction * numberOfParticles + 0.5));
        effectiveNumParticles += numParticlesInBin[k];
        if (numParticlesInBin[k] > numParticlesInBin[largestBin]) largestBin = k;
    }

    numParticlesInBin[largestBin] += (numberOfParticles - effectiveNumParticles);

    for (int k = 0; k < numberOfEnergyBins; ++ k) {
        gsl_ran_discrete_t *table
            = gsl_ran_discrete_preproc(numberOfSampleBins,
                                       &(distributionTable[numberOfSampleBins * k]));

        for (int i = 0; i < numParticlesInBin[k]; i++) {
            double xx[2];
            gsl_qrng_get(quasiRandGen, xx);
            tCoord = hi * (xx[1] + static_cast<int>(gsl_ran_discrete(randGen_m, table))
                           - binTotal / 2 + k * numberOfSampleBins) / (binTotal / 2);

            saveProcessor++;
            if (saveProcessor >= numNodes)
                saveProcessor = 0;

            if (scalable || myNode == saveProcessor) {
                tOrZDist_m.push_back(tCoord);
                pzDist_m.push_back(0.0);
            }
        }
        gsl_ran_discrete_free(table);
    }

    gsl_qrng_free(quasiRandGen);
    delete [] distributionTable;

}

void Distribution::generateBinomial(size_t numberOfParticles) {

    /*!
     *
     * \brief Following W. Johos for his report  <a href="https://intranet.psi.ch/pub/AUTHOR_WWW/ABE/TalksDE/TM-11-14.pdf"> TM-11-14 </a>
     *
     * For the \f$x,p_x\f$ phase space we have:
     * \f[
     *  \epsilon_x = \sigma_x \sigma_{p_x} \cos{( \arcsin{(\sigma_{12}) }) }
     *  \f]
     *
     * \f{eqnarray*}{
     * \beta_x &=& \frac{\sigma_x^2}{\epsilon_x} \\
     * \gamma_x &=& \frac{\sigma_{p_x}^2}{\epsilon_x} \\
     * \alpha_x &=& -\sigma_{12} \sqrt{(\beta_x \gamma_x)} \\
     * \f}
     *
     * This holds similar for the other dimensions.
     */


    // Calculate emittance and Twiss parameters.
    Vector_t emittance;
    Vector_t alpha, beta, gamma;
    for (unsigned int index = 0; index < 3; index++) {
        emittance(index) = sigmaR_m[index] * sigmaP_m[index]
            * cos(asin(correlationMatrix_m(2 * index + 1, 2 * index)));

        if (std::abs(emittance(index)) > std::numeric_limits<double>::epsilon()) {
            beta(index)  = pow(sigmaR_m[index], 2.0) / emittance(index);
            gamma(index) = pow(sigmaP_m[index], 2.0) / emittance(index);
        } else {
            beta(index)  = sqrt(std::numeric_limits<double>::max());
            gamma(index) = sqrt(std::numeric_limits<double>::max());
        }
        alpha(index) = -correlationMatrix_m(2 * index + 1, 2 * index)
                        * sqrt(beta(index) * gamma(index));
    }

    Vector_t M, PM, L, PL, X, PX;
    Vector_t CHI, COSCHI, SINCHI(0.0);
    Vector_t AMI;
    Vektor<BinomialBehaviorSplitter*, 3> splitter;
    for (unsigned int index = 0; index < 3; index++) {
        // gamma(index) *= cutoffR_m[index];
        // beta(index)  *= cutoffP_m[index];
        COSCHI[index] =  sqrt(1.0 / (1.0 + pow(alpha(index), 2.0)));
        SINCHI[index] = -alpha(index) * COSCHI[index];
        CHI[index]    =  atan2(SINCHI[index], COSCHI[index]);
        AMI[index]    =  1.0 / mBinomial_m[index];
        M[index]      =  sqrt(emittance(index) * beta(index));
        PM[index]     =  sqrt(emittance(index) * gamma(index));
        L[index]      =  sqrt((mBinomial_m[index] + 1.0) / 2.0) * M[index];
        PL[index]     =  sqrt((mBinomial_m[index] + 1.0) / 2.0) * PM[index];

        if (mBinomial_m[index] < 10000) {
            X[index] =  L[index];
            PX[index] = PL[index];
            splitter[index] = new MDependentBehavior(mBinomial_m[index]);
        } else {
            X[index] =  M[index];
            PX[index] = PM[index];
            splitter[index] = new GaussianLikeBehavior();
        }
    }

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    double temp = 1.0 - std::pow(correlationMatrix_m(1, 0), 2.0);
    const double tempa = copysign(sqrt(std::abs(temp)), temp);
    const double l32 = (correlationMatrix_m(4, 1) -
                        correlationMatrix_m(1, 0) * correlationMatrix_m(4, 0)) / tempa;
    temp = 1 - std::pow(correlationMatrix_m(4, 0), 2.0) - l32 * l32;
    const double l33 = copysign(sqrt(std::abs(temp)), temp);

    const double l42 = (correlationMatrix_m(5, 1) -
                        correlationMatrix_m(1, 0) * correlationMatrix_m(5, 0)) / tempa;
    const double l43 = (correlationMatrix_m(5, 4) -
                        correlationMatrix_m(4, 0) * correlationMatrix_m(5, 0) - l42 * l32) / l33;
    temp = 1 - std::pow(correlationMatrix_m(5, 0), 2.0) - l42 * l42 - l43 * l43;
    const double l44 = copysign(sqrt(std::abs(temp)), temp);

    Vector_t x = Vector_t(0.0);
    Vector_t p = Vector_t(0.0);
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double A = 0.0;
        double AL = 0.0;
        double Ux = 0.0, U = 0.0;
        double Vx = 0.0, V = 0.0;

        A = splitter[0]->get(gsl_rng_uniform(randGen_m));
        AL = Physics::two_pi * gsl_rng_uniform(randGen_m);
        Ux = A * cos(AL);
        Vx = A * sin(AL);
        x[0] = X[0] * Ux;
        p[0] = PX[0] * (Ux * SINCHI[0] + Vx * COSCHI[0]);

        A = splitter[1]->get(gsl_rng_uniform(randGen_m));
        AL = Physics::two_pi * gsl_rng_uniform(randGen_m);
        U = A * cos(AL);
        V = A * sin(AL);
        x[1] = X[1] * U;
        p[1] = PX[1] * (U * SINCHI[1] + V * COSCHI[1]);

        A = splitter[2]->get(gsl_rng_uniform(randGen_m));
        AL = Physics::two_pi * gsl_rng_uniform(randGen_m);
        U = A * cos(AL);
        V = A * sin(AL);
        x[2] = X[2] *  (Ux * correlationMatrix_m(4, 0) + Vx * l32 + U * l33);
        p[2] = PX[2] * (Ux * correlationMatrix_m(5, 0) + Vx * l42 + U * l43 + V * l44);

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x[0]);
            pxDist_m.push_back(p[0]);
            yDist_m.push_back(x[1]);
            pyDist_m.push_back(p[1]);
            tOrZDist_m.push_back(x[2]);
            pzDist_m.push_back(avrgpz_m + p[2]);
        }
    }

    for (unsigned int index = 0; index < 3; index++) {
        delete splitter[index];
    }
}

void Distribution::generateFlattopLaserProfile(size_t numberOfParticles) {

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        laserProfile_m->getXY(x, y);

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x * sigmaR_m[0]);
            pxDist_m.push_back(px);
            yDist_m.push_back(y * sigmaR_m[1]);
            pyDist_m.push_back(py);
        }
    }

    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        generateAstraFlattopT(numberOfParticles);
    else
        generateLongFlattopT(numberOfParticles);

}

void Distribution::generateFlattopT(size_t numberOfParticles) {

    gsl_qrng *quasiRandGen2D = selectRandomGenerator(Options::rngtype,2);

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        bool allow = false;
        while (!allow) {

            if (quasiRandGen2D != NULL) {
                double randNums[2];
                gsl_qrng_get(quasiRandGen2D, randNums);
                x = -1.0 + 2.0 * randNums[0];
                y = -1.0 + 2.0 * randNums[1];
            } else {
                x = -1.0 + 2.0 * gsl_rng_uniform(randGen_m);
                y = -1.0 + 2.0 * gsl_rng_uniform(randGen_m);
            }

            allow = (pow(x, 2.0) + pow(y, 2.0) <= 1.0);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
        }

    }

    gsl_qrng_free(quasiRandGen2D);

    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        generateAstraFlattopT(numberOfParticles);
    else
        generateLongFlattopT(numberOfParticles);

}

void Distribution::generateFlattopZ(size_t numberOfParticles) {

    gsl_qrng *quasiRandGen1D = selectRandomGenerator(Options::rngtype,1);
    gsl_qrng *quasiRandGen2D = selectRandomGenerator(Options::rngtype,2);

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;
        double z = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {

            if (quasiRandGen2D != NULL) {
                double randNums[2];
                gsl_qrng_get(quasiRandGen2D, randNums);
                x = -1.0 + 2.0 * randNums[0];
                y = -1.0 + 2.0 * randNums[1];
            } else {
                x = -1.0 + 2.0 * gsl_rng_uniform(randGen_m);
                y = -1.0 + 2.0 * gsl_rng_uniform(randGen_m);
            }

            allow = (pow(x, 2.0) + pow(y, 2.0) <= 1.0);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];

        if (quasiRandGen1D != NULL)
            gsl_qrng_get(quasiRandGen1D, &z);
        else
            z = gsl_rng_uniform(randGen_m);

        z = (z - 0.5) * sigmaR_m[2];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
            tOrZDist_m.push_back(z);
            pzDist_m.push_back(avrgpz_m + pz);
        }
    }

    gsl_qrng_free(quasiRandGen1D);
    gsl_qrng_free(quasiRandGen2D);
}

void Distribution::generateGaussZ(size_t numberOfParticles) {

    gsl_matrix * corMat  = gsl_matrix_alloc (6, 6);
    gsl_vector * rx = gsl_vector_alloc(6);
    gsl_vector * ry = gsl_vector_alloc(6);

    for (unsigned int i = 0; i < 6; ++ i) {
        gsl_matrix_set(corMat, i, i, correlationMatrix_m(i, i));
        for (unsigned int j = 0; j < i; ++ j) {
            gsl_matrix_set(corMat, i, j, correlationMatrix_m(i, j));
            gsl_matrix_set(corMat, j, i, correlationMatrix_m(i, j));
        }
    }

#define DISTDBG1
#ifdef DISTDBG1
    *gmsg << "* m before gsl_linalg_cholesky_decomp" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (j==0)
                *gmsg << "* r(" << std::setprecision(1) << i << ", : ) = "
                      << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
            else
                *gmsg << " & " << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
        }
        *gmsg << " \\\\" << endl;
    }
#endif
/*
    //Sets the GSL error handler off, exception will be handled internaly with a renormalization method
    gsl_set_error_handler_off();
*/
    int errcode = gsl_linalg_cholesky_decomp(corMat);
/*
    double rn = 1e-12;

    while (errcode == GSL_EDOM) {

        // Resets the correlation matrix
        for (unsigned int i = 0; i < 6; ++ i) {
            gsl_matrix_set(corMat, i, i, correlationMatrix_m(i, i));
            for (unsigned int j = 0; j < i; ++ j) {
                gsl_matrix_set(corMat, i, j, correlationMatrix_m(i, j));
                gsl_matrix_set(corMat, j, i, correlationMatrix_m(i, j));
            }
        }
        // Applying a renormalization method corMat = corMat + rn*Unitymatrix
        // This is the renormalization
        for(int i = 0; i < 6; i++){
            double corMati = gsl_matrix_get(corMat, i , i);
            gsl_matrix_set(corMat, i, i, corMati + rn);
        }
        //Just to be sure that the renormalization worked
        errcode = gsl_linalg_cholesky_decomp(corMat);
        if(errcode != GSL_EDOM) *gmsg << "* WARNING: Correlation matrix had to be renormalized"<<endl;
        else rn *= 10;
    }
    //Sets again the standard GSL error handler on
    gsl_set_error_handler(NULL);
*/
    //Just to be sure
    if (errcode == GSL_EDOM) {
        throw OpalException("Distribution::GenerateGaussZ",
                            "gsl_linalg_cholesky_decomp failed");
    }
    // so make the upper part zero.
    for (int i = 0; i < 5; ++ i) {
        for (int j = i+1; j < 6 ; ++ j) {
            gsl_matrix_set (corMat, i, j, 0.0);
        }
    }
#define DISTDBG2
#ifdef DISTDBG2
    *gmsg << "* m after gsl_linalg_cholesky_decomp" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (j==0)
                *gmsg << "* r(" << std::setprecision(1) << i << ", : ) = "
                      << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
            else
                *gmsg << " & " << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
        }
        *gmsg << " \\\\" << endl;
    }
#endif

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {
        bool allow = false;

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;
        double z = 0.0;
        double pz = 0.0;

        while (!allow) {
            gsl_vector_set(rx, 0, gsl_ran_gaussian(randGen_m, 1.0));
            gsl_vector_set(rx, 1, gsl_ran_gaussian(randGen_m, 1.0));
            gsl_vector_set(rx, 2, gsl_ran_gaussian(randGen_m, 1.0));
            gsl_vector_set(rx, 3, gsl_ran_gaussian(randGen_m, 1.0));
            gsl_vector_set(rx, 4, gsl_ran_gaussian(randGen_m, 1.0));
            gsl_vector_set(rx, 5, gsl_ran_gaussian(randGen_m, 1.0));

            gsl_blas_dgemv(CblasNoTrans, 1.0, corMat, rx, 0.0, ry);

            x = gsl_vector_get(ry, 0);
            px = gsl_vector_get(ry, 1);
            y = gsl_vector_get(ry, 2);
            py = gsl_vector_get(ry, 3);
            z = gsl_vector_get(ry, 4);
            pz = gsl_vector_get(ry, 5);

            bool xAndYOk = false;
            if (cutoffR_m[0] < SMALLESTCUTOFF && cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = true;
            else if (cutoffR_m[0] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(y) <= cutoffR_m[1]);
            else if (cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(x) <= cutoffR_m[0]);
            else
                xAndYOk = (pow(x / cutoffR_m[0], 2.0) + pow(y / cutoffR_m[1], 2.0) <= 1.0);

            bool pxAndPyOk = false;
            if (cutoffP_m[0] < SMALLESTCUTOFF && cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = true;
            else if (cutoffP_m[0] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(py) <= cutoffP_m[1]);
            else if (cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(px) <= cutoffP_m[0]);
            else
                pxAndPyOk = (pow(px / cutoffP_m[0], 2.0) + pow(py / cutoffP_m[1], 2.0) <= 1.0);

            allow = (xAndYOk && pxAndPyOk && std::abs(z) < cutoffR_m[2] && std::abs(pz) < cutoffP_m[2]);
        }

        x *= sigmaR_m[0];
        y *= sigmaR_m[1];
        z *= sigmaR_m[2];
        px *= sigmaP_m[0];
        py *= sigmaP_m[1];
        pz *= sigmaP_m[2];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
            tOrZDist_m.push_back(z);
            pzDist_m.push_back(avrgpz_m + pz);
        }
    }

    gsl_vector_free(rx);
    gsl_vector_free(ry);
    gsl_matrix_free(corMat);
}

void Distribution::generateLongFlattopT(size_t numberOfParticles) {

    double flattopTime = tPulseLengthFWHM_m
        - sqrt(2.0 * log(2.0)) * (sigmaTRise_m + sigmaTFall_m);

    if (flattopTime < 0.0)
        flattopTime = 0.0;

    double normalizedFlankArea = 0.5 * sqrt(2.0 * Physics::pi) * gsl_sf_erf(cutoffR_m[2] / sqrt(2.0));
    double distArea = flattopTime
        + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea;

    // Find number of particles in rise, fall and flat top.
    size_t numRise = numberOfParticles * sigmaTRise_m * normalizedFlankArea / distArea;
    size_t numFall = numberOfParticles * sigmaTFall_m * normalizedFlankArea / distArea;
    size_t numFlat = numberOfParticles - numRise - numFall;

    // Generate particles in tail.
    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    for (size_t partIndex = 0; partIndex < numFall; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {
            t = gsl_ran_gaussian_tail(randGen_m, 0, sigmaTFall_m);
            if (t <= sigmaTRise_m * cutoffR_m[2]) {
                t = -t + sigmaTFall_m * cutoffR_m[2];
                allow = true;
            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    /*
     * Generate particles in flat top. The flat top can also have sinusoidal
     * modulations.
     */
    double modulationAmp = Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCAMPLITUDE]) / 100.0;
    if (modulationAmp > 1.0)
        modulationAmp = 1.0;
    double numModulationPeriods
        = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCPERIODS]));
    double modulationPeriod = 0.0;
    if (numModulationPeriods != 0.0)
        modulationPeriod = flattopTime / numModulationPeriods;

    gsl_qrng *quasiRandGen1D = selectRandomGenerator(Options::rngtype,1);
    gsl_qrng *quasiRandGen2D = selectRandomGenerator(Options::rngtype,2);

    for (size_t partIndex = 0; partIndex < numFlat; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        if (modulationAmp == 0.0 || numModulationPeriods == 0.0) {

            if (quasiRandGen1D != NULL)
                gsl_qrng_get(quasiRandGen1D, &t);
            else
                t = gsl_rng_uniform(randGen_m);

            t = flattopTime * t + sigmaTFall_m * cutoffR_m[2];

        } else {

            double funcAmp = 0.0;
            bool allow = false;
            while (!allow) {
                if (quasiRandGen2D != NULL) {
                    double randNums[2] = {0.0, 0.0};
                    gsl_qrng_get(quasiRandGen2D, randNums);
                    t = randNums[0] * flattopTime;
                    funcAmp = randNums[1];
                } else {
                    t = gsl_rng_uniform(randGen_m)*flattopTime;
                    funcAmp = gsl_rng_uniform(randGen_m);
                }

                double funcValue = (1.0 + modulationAmp
                                    * sin(Physics::two_pi * t / modulationPeriod))
                    / (1.0 + std::abs(modulationAmp));

                allow = (funcAmp <= funcValue);

            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    // Generate particles in rise.
    for (size_t partIndex = 0; partIndex < numRise; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {
            t = gsl_ran_gaussian_tail(randGen_m, 0, sigmaTRise_m);
            if (t <= sigmaTRise_m * cutoffR_m[2]) {
                t += sigmaTFall_m * cutoffR_m[2] + flattopTime;
                allow = true;
            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    gsl_qrng_free(quasiRandGen1D);
    gsl_qrng_free(quasiRandGen2D);
}

void Distribution::generateTransverseGauss(size_t numberOfParticles) {

    // Generate coordinates.
    gsl_matrix * corMat  = gsl_matrix_alloc (4, 4);
    gsl_vector * rx = gsl_vector_alloc(4);
    gsl_vector * ry = gsl_vector_alloc(4);

    for (unsigned int i = 0; i < 4; ++ i) {
        gsl_matrix_set(corMat, i, i, correlationMatrix_m(i, i));
        for (unsigned int j = 0; j < i; ++ j) {
            gsl_matrix_set(corMat, i, j, correlationMatrix_m(i, j));
            gsl_matrix_set(corMat, j, i, correlationMatrix_m(i, j));
        }
    }

    int errcode = gsl_linalg_cholesky_decomp(corMat);

    if (errcode == GSL_EDOM) {
        throw OpalException("Distribution::GenerateTransverseGauss",
                            "gsl_linalg_cholesky_decomp failed");
    }

    for (int i = 0; i < 3; ++ i) {
        for (int j = i+1; j < 4 ; ++ j) {
            gsl_matrix_set (corMat, i, j, 0.0);
        }
    }

    int saveProcessor = -1;
    const int myNode = Ippl::myNode();
    const int numNodes = Ippl::getNodes();
    const bool scalable = Attributes::getBool(itsAttr[Attrib::Distribution::SCALABLE]);

    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        bool allow = false;
        while (!allow) {

            gsl_vector_set(rx, 0, gsl_ran_gaussian (randGen_m,1.0));
            gsl_vector_set(rx, 1, gsl_ran_gaussian (randGen_m,1.0));
            gsl_vector_set(rx, 2, gsl_ran_gaussian (randGen_m,1.0));
            gsl_vector_set(rx, 3, gsl_ran_gaussian (randGen_m,1.0));

            gsl_blas_dgemv(CblasNoTrans, 1.0, corMat, rx, 0.0, ry);
            x = gsl_vector_get(ry, 0);
            px = gsl_vector_get(ry, 1);
            y = gsl_vector_get(ry, 2);
            py = gsl_vector_get(ry, 3);

            bool xAndYOk = false;
            if (cutoffR_m[0] < SMALLESTCUTOFF && cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = true;
            else if (cutoffR_m[0] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(y) <= cutoffR_m[1]);
            else if (cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(x) <= cutoffR_m[0]);
            else
                xAndYOk = (pow(x / cutoffR_m[0], 2.0) + pow(y / cutoffR_m[1], 2.0) <= 1.0);

            bool pxAndPyOk = false;
            if (cutoffP_m[0] < SMALLESTCUTOFF && cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = true;
            else if (cutoffP_m[0] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(py) <= cutoffP_m[1]);
            else if (cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(px) <= cutoffP_m[0]);
            else
                pxAndPyOk = (pow(px / cutoffP_m[0], 2.0) + pow(py / cutoffP_m[1], 2.0) <= 1.0);

            allow = (xAndYOk && pxAndPyOk);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];
        px *= sigmaP_m[0];
        py *= sigmaP_m[1];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= numNodes)
            saveProcessor = 0;

        if (scalable || myNode == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
        }
    }

    gsl_matrix_free(corMat);
    gsl_vector_free(rx);
    gsl_vector_free(ry);
}

void Distribution::initializeBeam(PartBunchBase<double, 3> *beam) {

    /*
     * Set emission time, the current beam bunch number and
     * set the beam energy bin structure, if it has one.
     */
    beam->setTEmission(tEmission_m);
    beam->setNumBunch(1);
    if (numberOfEnergyBins_m > 0)
        beam->setEnergyBins(numberOfEnergyBins_m);
}

void Distribution::injectBeam(PartBunchBase<double, 3> *beam) {

    writeOutFileInjection();

    std::vector<double> id1 = Attributes::getRealArray(itsAttr[Attrib::Distribution::ID1]);
    std::vector<double> id2 = Attributes::getRealArray(itsAttr[Attrib::Distribution::ID2]);

    bool hasID1 = (id1.size() != 0);
    bool hasID2 = (id2.size() != 0);

    if (hasID1 || hasID2)
        *gmsg << "* Use special ID1 or ID2 particle in distribution" << endl;


    size_t numberOfParticles = tOrZDist_m.size();
    beam->create(numberOfParticles);
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        beam->R[partIndex] = Vector_t(xDist_m.at(partIndex),
                                     yDist_m.at(partIndex),
                                     tOrZDist_m.at(partIndex));

        beam->P[partIndex] = Vector_t(pxDist_m.at(partIndex),
                                      pyDist_m.at(partIndex),
                                      pzDist_m.at(partIndex));

        beam->Q[partIndex] = beam->getChargePerParticle();
        beam->Ef[partIndex] = Vector_t(0.0);
        beam->Bf[partIndex] = Vector_t(0.0);
        beam->PType[partIndex] = ParticleType::REGULAR;
        beam->TriID[partIndex] = 0;

        if (numberOfEnergyBins_m > 0) {
            size_t binNumber = findEBin(tOrZDist_m.at(partIndex));
            beam->Bin[partIndex] = binNumber;
            beam->iterateEmittedBin(binNumber);
        } else
            beam->Bin[partIndex] = 0;

        if (hasID1) {
            if (beam->ID[partIndex] == 1) {
                beam->R[partIndex] = Vector_t(id1[0],id1[1],id1[2]);
                beam->P[partIndex] = Vector_t(id1[3],id1[4],id1[5]);
            }
        }

        if (hasID2) {
            if (beam->ID[partIndex] == 2) {
                beam->R[partIndex] = Vector_t(id2[0],id2[1],id2[2]);
                beam->P[partIndex] = Vector_t(id2[3],id2[4],id2[5]);
            }
        }
    }

    xDist_m.clear();
    pxDist_m.clear();
    yDist_m.clear();
    pyDist_m.clear();
    tOrZDist_m.clear();
    pzDist_m.clear();

    beam->boundp();
    beam->calcEMean();
}

bool Distribution::getIfDistEmitting() {
    return emitting_m;
}

int Distribution::getLastEmittedEnergyBin() {
    return currentEnergyBin_m;
}

double Distribution::getMaxTOrZ() {

    std::vector<double>::iterator longIt = tOrZDist_m.begin();
    double maxTOrZ = *longIt;
    for (++longIt; longIt != tOrZDist_m.end(); ++longIt) {
        if (maxTOrZ < *longIt)
            maxTOrZ = *longIt;
    }

    reduce(maxTOrZ, maxTOrZ, OpMaxAssign());

    return maxTOrZ;
}

double Distribution::getMinTOrZ() {

    std::vector<double>::iterator longIt = tOrZDist_m.begin();
    double minTOrZ = *longIt;
    for (++longIt; longIt != tOrZDist_m.end(); ++longIt) {
        if (minTOrZ > *longIt)
            minTOrZ = *longIt;
    }

    reduce(minTOrZ, minTOrZ, OpMinAssign());

    return minTOrZ;
}

size_t Distribution::getNumberOfEmissionSteps() {
    return static_cast<size_t> (numberOfEnergyBins_m * numberOfSampleBins_m);
}

int Distribution::getNumberOfEnergyBins() {
    return numberOfEnergyBins_m;
}

double Distribution::getEmissionDeltaT() {
    return tBin_m / numberOfSampleBins_m;
}

double Distribution::getEnergyBinDeltaT() {
    return tBin_m;
}

double Distribution::getWeight() {
    return std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::WEIGHT]));
}

std::vector<double>& Distribution::getXDist() {
    return xDist_m;
}

std::vector<double>& Distribution::getBGxDist() {
    return pxDist_m;
}

std::vector<double>& Distribution::getYDist() {
    return yDist_m;
}

std::vector<double>& Distribution::getBGyDist() {
    return pyDist_m;
}

std::vector<double>& Distribution::getTOrZDist() {
    return tOrZDist_m;
}

std::vector<double>& Distribution::getBGzDist() {
    return pzDist_m;
}

void Distribution::printDist(Inform &os, size_t numberOfParticles) const {

    if (numberOfParticles > 0) {
        size_t np = numberOfParticles * (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE)? 2: 1);
        reduce(np, np, OpAddAssign());
        os << "* Number of particles: "
           << np
           << endl
           << "* " << endl;
    }

    switch (distrTypeT_m) {

    case DistrTypeT::FROMFILE:
        printDistFromFile(os);
        break;
    case DistrTypeT::GAUSS:
        printDistGauss(os);
        break;
    case DistrTypeT::BINOMIAL:
        printDistBinomial(os);
        break;
    case DistrTypeT::SURFACEEMISSION:
        printDistSurfEmission(os);
        break;
    case DistrTypeT::SURFACERANDCREATE:
        printDistSurfAndCreate(os);
        break;
    case DistrTypeT::FLATTOP:
    case DistrTypeT::GUNGAUSSFLATTOPTH:
    case DistrTypeT::ASTRAFLATTOPTH:
        printDistFlattop(os);
        break;
    case DistrTypeT::MATCHEDGAUSS:
        printDistMatchedGauss(os);
        break;
    default:
        INFOMSG("Distribution unknown." << endl;);
        break;
    }

}

void Distribution::printDistBinomial(Inform &os) const {

    os << "* Distribution type: BINOMIAL" << endl;
    os << "* " << endl;
    os << "* SIGMAX   = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY   = " << sigmaR_m[1] << " [m]" << endl;
    if (emitting_m)
        os << "* SIGMAT   = " << sigmaR_m[2] << " [sec]" << endl;
    else
        os << "* SIGMAZ   = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX  = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY  = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ  = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
    os << "* MX       = " << mBinomial_m[0] << endl;
    os << "* MY       = " << mBinomial_m[1] << endl;
    if (emitting_m)
        os << "* MT       = " << mBinomial_m[2] << endl;
    else
        os << "* MZ       = " << mBinomial_m[2] << endl;
    os << "* CORRX    = " << correlationMatrix_m(1, 0) << endl;
    os << "* CORRY    = " << correlationMatrix_m(3, 2) << endl;
    os << "* CORRZ    = " << correlationMatrix_m(5, 4) << endl;
    os << "* R61      = " << correlationMatrix_m(5, 0) << endl;
    os << "* R62      = " << correlationMatrix_m(5, 1) << endl;
    os << "* R51      = " << correlationMatrix_m(4, 0) << endl;
    os << "* R52      = " << correlationMatrix_m(4, 1) << endl;
}

void Distribution::printDistFlattop(Inform &os) const {

    switch (distrTypeT_m) {

    case DistrTypeT::ASTRAFLATTOPTH:
        os << "* Distribution type: ASTRAFLATTOPTH" << endl;
        break;

    case DistrTypeT::GUNGAUSSFLATTOPTH:
        os << "* Distribution type: GUNGAUSSFLATTOPTH" << endl;
        break;

    default:
        os << "* Distribution type: FLATTOP" << endl;
        break;

    }
    os << "* " << endl;

    if (laserProfile_m != NULL) {

        os << "* Transverse profile determined by laser image: " << endl
           << endl
           << "* Laser profile: " << laserProfileFileName_m << endl
           << "* Laser image:   " << laserImageName_m << endl
           << "* Intensity cut: " << laserIntensityCut_m << endl;

    } else {

        os << "* SIGMAX                        = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY                        = " << sigmaR_m[1] << " [m]" << endl;

    }

    if (emitting_m) {

        if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH) {

            os << "* Time Rise                     = " << tRise_m
               << " [sec]" << endl;
            os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
               << " [sec]" << endl;

        } else {
            os << "* Sigma Time Rise               = " << sigmaTRise_m
               << " [sec]" << endl;
            os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
               << " [sec]" << endl;
            os << "* Sigma Time Fall               = " << sigmaTFall_m
               << " [sec]" << endl;
            os << "* Longitudinal cutoff           = " << cutoffR_m[2]
               << " [units of Sigma Time]" << endl;
            os << "* Flat top modulation amplitude = "
               << Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCAMPLITUDE])
               << " [Percent of distribution amplitude]" << endl;
            os << "* Flat top modulation periods   = "
               << std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCPERIODS]))
               << endl;
        }

    } else
        os << "* SIGMAZ                        = " << sigmaR_m[2]
           << " [m]" << endl;
}

void Distribution::printDistFromFile(Inform &os) const {
    os << "* Distribution type: FROMFILE" << endl;
    os << "* " << endl;
    os << "* Input file:        "
       << Attributes::getString(itsAttr[Attrib::Distribution::FNAME]) << endl;
}


void Distribution::printDistMatchedGauss(Inform &os) const {
    os << "* Distribution type: MATCHEDGAUSS" << endl;
    os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;
    os << "* SIGMAZ     = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX    = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY    = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ    = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
    os << "* AVRGPZ     = " << avrgpz_m <<    " [Beta Gamma]" << endl;

    os << "* CORRX      = " << correlationMatrix_m(1, 0) << endl;
    os << "* CORRY      = " << correlationMatrix_m(3, 2) << endl;
    os << "* CORRZ      = " << correlationMatrix_m(5, 4) << endl;
    os << "* R61        = " << correlationMatrix_m(5, 0) << endl;
    os << "* R62        = " << correlationMatrix_m(5, 1) << endl;
    os << "* R51        = " << correlationMatrix_m(4, 0) << endl;
    os << "* R52        = " << correlationMatrix_m(4, 1) << endl;
    os << "* CUTOFFX    = " << cutoffR_m[0] << " [units of SIGMAX]" << endl;
    os << "* CUTOFFY    = " << cutoffR_m[1] << " [units of SIGMAY]" << endl;
    os << "* CUTOFFLONG = " << cutoffR_m[2] << " [units of SIGMAZ]" << endl;
    os << "* CUTOFFPX   = " << cutoffP_m[0] << " [units of SIGMAPX]" << endl;
    os << "* CUTOFFPY   = " << cutoffP_m[1] << " [units of SIGMAPY]" << endl;
    os << "* CUTOFFPZ   = " << cutoffP_m[2] << " [units of SIGMAPY]" << endl;
}

void Distribution::printDistGauss(Inform &os) const {
    os << "* Distribution type: GAUSS" << endl;
    os << "* " << endl;
    if (emitting_m) {
        os << "* Sigma Time Rise               = " << sigmaTRise_m
           << " [sec]" << endl;
        os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
           << " [sec]" << endl;
        os << "* Sigma Time Fall               = " << sigmaTFall_m
           << " [sec]" << endl;
        os << "* Longitudinal cutoff           = " << cutoffR_m[2]
           << " [units of Sigma Time]" << endl;
        os << "* Flat top modulation amplitude = "
           << Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCAMPLITUDE])
           << " [Percent of distribution amplitude]" << endl;
        os << "* Flat top modulation periods   = "
           << std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::FTOSCPERIODS]))
           << endl;
        os << "* SIGMAX                        = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY                        = " << sigmaR_m[1] << " [m]" << endl;
        os << "* SIGMAPX                       = " << sigmaP_m[0]
           << " [Beta Gamma]" << endl;
        os << "* SIGMAPY                       = " << sigmaP_m[1]
           << " [Beta Gamma]" << endl;
        os << "* CORRX                         = " << correlationMatrix_m(1, 0) << endl;
        os << "* CORRY                         = " << correlationMatrix_m(3, 2) << endl;
        os << "* CUTOFFX                       = " << cutoffR_m[0]
           << " [units of SIGMAX]" << endl;
        os << "* CUTOFFY                       = " << cutoffR_m[1]
           << " [units of SIGMAY]" << endl;
        os << "* CUTOFFPX                      = " << cutoffP_m[0]
           << " [units of SIGMAPX]" << endl;
        os << "* CUTOFFPY                      = " << cutoffP_m[1]
           << " [units of SIGMAPY]" << endl;
    } else {
        os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;
        os << "* SIGMAZ     = " << sigmaR_m[2] << " [m]" << endl;
        os << "* SIGMAPX    = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
        os << "* SIGMAPY    = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
        os << "* SIGMAPZ    = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
        os << "* AVRGPZ     = " << avrgpz_m <<    " [Beta Gamma]" << endl;

        os << "* CORRX      = " << correlationMatrix_m(1, 0) << endl;
        os << "* CORRY      = " << correlationMatrix_m(3, 2) << endl;
        os << "* CORRZ      = " << correlationMatrix_m(5, 4) << endl;
        os << "* R61        = " << correlationMatrix_m(5, 0) << endl;
        os << "* R62        = " << correlationMatrix_m(5, 1) << endl;
        os << "* R51        = " << correlationMatrix_m(4, 0) << endl;
        os << "* R52        = " << correlationMatrix_m(4, 1) << endl;
        os << "* CUTOFFX    = " << cutoffR_m[0] << " [units of SIGMAX]" << endl;
        os << "* CUTOFFY    = " << cutoffR_m[1] << " [units of SIGMAY]" << endl;
        os << "* CUTOFFLONG = " << cutoffR_m[2] << " [units of SIGMAZ]" << endl;
        os << "* CUTOFFPX   = " << cutoffP_m[0] << " [units of SIGMAPX]" << endl;
        os << "* CUTOFFPY   = " << cutoffP_m[1] << " [units of SIGMAPY]" << endl;
        os << "* CUTOFFPZ   = " << cutoffP_m[2] << " [units of SIGMAPY]" << endl;
    }
}

void Distribution::printDistSurfEmission(Inform &os) const {

    os << "* Distribution type: SURFACEEMISION" << endl;
    os << "* " << endl;
    os << "* * Number of electrons for surface emission  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::NPDARKCUR]) << endl;
    os << "* * Initialized electrons inward margin for surface emission  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::INWARDMARGIN]) << endl;
    os << "* * E field threshold (MV), only in position r with E(r)>EINITHR that "
       << "particles will be initialized   "
       << Attributes::getReal(itsAttr[Attrib::Distribution::EINITHR]) << endl;
    os << "* * Field enhancement for surface emission  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::FNBETA]) << endl;
    os << "* * Maximum number of electrons emitted from a single triangle for "
       << "Fowler-Nordheim emission  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::FNMAXEMI]) << endl;
    os << "* * Field Threshold for Fowler-Nordheim emission (MV/m)  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNFIELDTHR]) << endl;
    os << "* * Empirical constant A for Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNA]) << endl;
    os << "* * Empirical constant B for Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNB]) << endl;
    os << "* * Constant for image charge effect parameter y(E) in Fowler-Nordheim "
       << "emission model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNY]) << endl;
    os << "* * Zero order constant for image charge effect parameter v(y) in "
       << "Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNVYZERO]) << endl;
    os << "* * Second order constant for image charge effect parameter v(y) in "
       << "Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::FNVYSECOND]) << endl;
    os << "* * Select secondary model type(0:no secondary emission; 1:Furman-Pivi; "
       << "2 or larger: Vaughan's model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::SECONDARYFLAG]) << endl;
    os << "* * Secondary emission mode type(true: emit n true secondaries; false: "
       << "emit one particle with n times charge  "
       << Attributes::getBool(itsAttr[Attrib::Distribution::NEMISSIONMODE]) << endl;
    os << "* * Sey_0 in Vaughan's model "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VSEYZERO]) << endl;
    os << "* * Energy related to sey_0 in Vaughan's model in eV  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VEZERO]) << endl;
    os << "* * Sey max in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VSEYMAX]) << endl;
    os << "* * Emax in Vaughan's model in eV  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VEMAX]) << endl;
    os << "* * Fitting parameter denotes the roughness of surface for impact "
       << "energy in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VKENERGY]) << endl;
    os << "* * Fitting parameter denotes the roughness of surface for impact angle "
       << "in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VKTHETA]) << endl;
    os << "* * Thermal velocity of Maxwellian distribution of secondaries in "
       << "Vaughan's model  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VVTHERMAL]) << endl;
    os << "* * VW denote the velocity scalar for Parallel plate benchmark  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::VW]) << endl;
    os << "* * Material type number of the cavity surface for Furman-Pivi's model, "
       << "0 for copper, 1 for stainless steel  "
       <<  Attributes::getReal(itsAttr[Attrib::Distribution::SURFMATERIAL]) << endl;
}

void Distribution::printDistSurfAndCreate(Inform &os) const {

    os << "* Distribution type: SURFACERANDCREATE" << endl;
    os << "* " << endl;
    os << "* * Number of electrons initialized on the surface as primaries  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::NPDARKCUR]) << endl;
    os << "* * Initialized electrons inward margin for surface emission  "
       << Attributes::getReal(itsAttr[Attrib::Distribution::INWARDMARGIN]) << endl;
    os << "* * E field threshold (MV), only in position r with E(r)>EINITHR that "
       << "particles will be initialized   "
       << Attributes::getReal(itsAttr[Attrib::Distribution::EINITHR]) << endl;
}

void Distribution::printEmissionModel(Inform &os) const {

    os << "* ------------- THERMAL EMITTANCE MODEL --------------------------------------------" << endl;

    switch (emissionModel_m) {

    case EmissionModelT::NONE:
        printEmissionModelNone(os);
        break;
    case EmissionModelT::ASTRA:
        printEmissionModelAstra(os);
        break;
    case EmissionModelT::NONEQUIL:
        printEmissionModelNonEquil(os);
        break;
    default:
        break;
    }

    os << "* ----------------------------------------------------------------------------------" << endl;

}

void Distribution::printEmissionModelAstra(Inform &os) const {
    os << "*  THERMAL EMITTANCE in ASTRA MODE" << endl;
    os << "*  Kinetic energy (thermal emittance) = "
       << std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]))
       << " [eV]  " << endl;
}

void Distribution::printEmissionModelNone(Inform &os) const {
    os << "*  THERMAL EMITTANCE in NONE MODE" << endl;
    os << "*  Kinetic energy added to longitudinal momentum = "
       << std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]))
       << " [eV]  " << endl;
}

void Distribution::printEmissionModelNonEquil(Inform &os) const {
    os << "*  THERMAL EMITTANCE in NONEQUIL MODE" << endl;
    os << "*  Cathode work function     = " << cathodeWorkFunc_m << " [eV]  " << endl;
    os << "*  Cathode Fermi energy      = " << cathodeFermiEnergy_m << " [eV]  " << endl;
    os << "*  Cathode temperature       = " << cathodeTemp_m / Physics::kB << " [K]  " << endl;
    os << "*  Photocathode laser energy = " << laserEnergy_m << " [eV]  " << endl;
}

void Distribution::printEnergyBins(Inform &os) const {

    os << "* " << endl;
    for (int binIndex = numberOfEnergyBins_m - 1; binIndex >=0; binIndex--) {
        size_t sum = 0;
        for (int sampleIndex = 0; sampleIndex < numberOfSampleBins_m; sampleIndex++)
            sum += gsl_histogram_get(energyBinHist_m,
                                     binIndex * numberOfSampleBins_m + sampleIndex);

        os << "* Energy Bin # " << numberOfEnergyBins_m - binIndex
           << " contains " << sum << " particles" << endl;
    }
    os << "* " << endl;

}

bool Distribution::Rebin() {
    /*
     * Allow a rebin (numberOfEnergyBins_m = 0) if all particles are emitted or
     * injected.
     */
    if (!emitting_m) {
        numberOfEnergyBins_m = 0;
        return true;
    } else {
        return false;
    }
}

void Distribution::reflectDistribution(size_t &numberOfParticles) {

    if (!Options::cZero || (distrTypeT_m == DistrTypeT::FROMFILE))
        return;

    size_t currentNumPart = tOrZDist_m.size();
    for (size_t partIndex = 0; partIndex < currentNumPart; partIndex++) {
        xDist_m.push_back(-xDist_m.at(partIndex));
        pxDist_m.push_back(-pxDist_m.at(partIndex));
        yDist_m.push_back(-yDist_m.at(partIndex));
        pyDist_m.push_back(-pyDist_m.at(partIndex));
        tOrZDist_m.push_back(tOrZDist_m.at(partIndex));
        pzDist_m.push_back(pzDist_m.at(partIndex));
    }
    numberOfParticles = tOrZDist_m.size();
    reduce(numberOfParticles, numberOfParticles, OpAddAssign());

    // if numberOfParticles is odd then delete last particle
    if (Ippl::myNode() == 0 &&
        (numberOfParticles + 1) / 2 != numberOfParticles / 2) {
        xDist_m.pop_back();
        pxDist_m.pop_back();
        yDist_m.pop_back();
        pyDist_m.pop_back();
        tOrZDist_m.pop_back();
        pzDist_m.pop_back();
    }
}

void Distribution::scaleDistCoordinates() {
    // at this point the distributions of an array of distributions are still separated

    const double xmult    = Attributes::getReal(itsAttr[Attrib::Distribution::XMULT]);
    const double pxmult   = Attributes::getReal(itsAttr[Attrib::Distribution::PXMULT]);
    const double ymult    = Attributes::getReal(itsAttr[Attrib::Distribution::YMULT]);
    const double pymult   = Attributes::getReal(itsAttr[Attrib::Distribution::PYMULT]);
    const double longmult = (emitting_m ?
                             Attributes::getReal(itsAttr[Attrib::Distribution::TMULT]) :
                             Attributes::getReal(itsAttr[Attrib::Distribution::ZMULT]));
    const double pzmult   = Attributes::getReal(itsAttr[Attrib::Distribution::PZMULT]);

    for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); ++ particleIndex) {
        xDist_m.at(particleIndex)    *= xmult;
        pxDist_m.at(particleIndex)   *= pxmult;
        yDist_m.at(particleIndex)    *= ymult;
        pyDist_m.at(particleIndex)   *= pymult;
        tOrZDist_m.at(particleIndex) *= longmult;
        pzDist_m.at(particleIndex)   *= pzmult;
    }
}

gsl_qrng* Distribution::selectRandomGenerator(std::string,unsigned int dimension)
{
    gsl_qrng *quasiRandGen = nullptr;
    if (Options::rngtype != std::string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if (Options::rngtype == std::string("HALTON")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, dimension);
        } else if (Options::rngtype == std::string("SOBOL")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_sobol, dimension);
        } else if (Options::rngtype == std::string("NIEDERREITER")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_niederreiter_2, dimension);
        } else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, dimension);
        }
    }
    return quasiRandGen;
}

void Distribution::setAttributes() {
    itsAttr[Attrib::Distribution::TYPE]
        = Attributes::makeString("TYPE","Distribution type: "
                                 "FROMFILE, "
                                 "GAUSS, "
                                 "BINOMIAL, "
                                 "FLATTOP, "
                                 "SURFACEEMISSION, "
                                 "SURFACERANDCREATE, "
                                 "GUNGAUSSFLATTOPTH, "
                                 "ASTRAFLATTOPTH, "
                                 "GAUSSMATCHED");
    itsAttr[Attrib::Legacy::Distribution::DISTRIBUTION]
        = Attributes::makeString("DISTRIBUTION","This attribute isn't supported any more. Use TYPE instead");
    itsAttr[Attrib::Distribution::LINE]
        = Attributes::makeString("LINE", "Beamline that contains a cyclotron or ring ", "");
    itsAttr[Attrib::Distribution::EX]
        = Attributes::makeReal("EX", "Projected normalized emittance EX (m-rad), used to create matched distribution ", 1E-6);
    itsAttr[Attrib::Distribution::EY]
        = Attributes::makeReal("EY", "Projected normalized emittance EY (m-rad) used to create matched distribution ", 1E-6);
    itsAttr[Attrib::Distribution::ET]
        = Attributes::makeReal("ET", "Projected normalized emittance ET (m-rad) used to create matched distribution ", 1E-6);
    // itsAttr[Attrib::Distribution::E2]
    //     = Attributes::makeReal("E2", "If E2<Eb, we compute the tunes from the beams energy Eb to E2 with dE=0.25 MeV ", 0.0);
    itsAttr[Attrib::Distribution::RESIDUUM]
        = Attributes::makeReal("RESIDUUM", "Residuum for the closed orbit finder and sigma matrix generator ", 1e-8);
    itsAttr[Attrib::Distribution::MAXSTEPSCO]
        = Attributes::makeReal("MAXSTEPSCO", "Maximum steps used to find closed orbit ", 100);
    itsAttr[Attrib::Distribution::MAXSTEPSSI]
        = Attributes::makeReal("MAXSTEPSSI", "Maximum steps used to find matched distribution ",500);
    itsAttr[Attrib::Distribution::ORDERMAPS]
        = Attributes::makeReal("ORDERMAPS", "Order used in the field expansion ", 7);
    itsAttr[Attrib::Distribution::SECTOR]
        = Attributes::makeBool("SECTOR", "Match using single sector (true)"
                               " (false: using all sectors) (default: true)",
                               true);
    itsAttr[Attrib::Distribution::NSTEPS]
        = Attributes::makeReal("NSTEPS", "Number of integration steps of closed orbit finder (matched gauss)"
                               " (default: 720)", 720);

    itsAttr[Attrib::Distribution::RGUESS]
        = Attributes::makeReal("RGUESS", "Guess value of radius (m) for closed orbit finder ", -1);

    itsAttr[Attrib::Distribution::DENERGY]
        = Attributes::makeReal("DENERGY", "Energy step size for closed orbit finder [GeV]", 0.001);

    itsAttr[Attrib::Distribution::NSECTORS]
        = Attributes::makeReal("NSECTORS", "Number of sectors to average field, for closed orbit finder ", 1);

    itsAttr[Attrib::Distribution::FNAME]
        = Attributes::makeString("FNAME", "File for reading in 6D particle "
                                 "coordinates.", "");


    itsAttr[Attrib::Distribution::WRITETOFILE]
        = Attributes::makeBool("WRITETOFILE", "Write initial distribution to file.",
                               false);

    itsAttr[Attrib::Distribution::WEIGHT]
        = Attributes::makeReal("WEIGHT", "Weight of distribution when used in a "
                               "distribution list.", 1.0);

    itsAttr[Attrib::Distribution::INPUTMOUNITS]
        = Attributes::makeString("INPUTMOUNITS", "Tell OPAL what input units are for momentum."
                                 " Currently \"NONE\" or \"EV\".", "");

    // Attributes for beam emission.
    itsAttr[Attrib::Distribution::EMITTED]
        = Attributes::makeBool("EMITTED", "Emitted beam, from cathode, as opposed to "
                               "an injected beam.", false);
    itsAttr[Attrib::Distribution::EMISSIONSTEPS]
        = Attributes::makeReal("EMISSIONSTEPS", "Number of time steps to use during emission.",
                               1);
    itsAttr[Attrib::Distribution::EMISSIONMODEL]
        = Attributes::makeString("EMISSIONMODEL", "Model used to emit electrons from a "
                                 "photocathode.", "None");
    itsAttr[Attrib::Distribution::EKIN]
        = Attributes::makeReal("EKIN", "Kinetic energy used in ASTRA thermal emittance "
                               "model (eV). (Thermal energy added in with random "
                               "number generator.)", 1.0);
    itsAttr[Attrib::Distribution::ELASER]
        = Attributes::makeReal("ELASER", "Laser energy (eV) for photocathode "
                               "emission. (Default 255 nm light.)", 4.86);
    itsAttr[Attrib::Distribution::W]
        = Attributes::makeReal("W", "Workfunction of photocathode material (eV)."
                               " (Default atomically clean copper.)", 4.31);
    itsAttr[Attrib::Distribution::FE]
        = Attributes::makeReal("FE", "Fermi energy of photocathode material (eV)."
                               " (Default atomically clean copper.)", 7.0);
    itsAttr[Attrib::Distribution::CATHTEMP]
        = Attributes::makeReal("CATHTEMP", "Temperature of photocathode (K)."
                               " (Default 300 K.)", 300.0);
    itsAttr[Attrib::Distribution::NBIN]
        = Attributes::makeReal("NBIN", "Number of energy bins to use when doing a space "
                               "charge solve.", 0.0);

    // Parameters for shifting or scaling initial distribution.
    itsAttr[Attrib::Distribution::XMULT] = Attributes::makeReal("XMULT", "Multiplier for x dimension.", 1.0);
    itsAttr[Attrib::Distribution::YMULT] = Attributes::makeReal("YMULT", "Multiplier for y dimension.", 1.0);
    itsAttr[Attrib::Distribution::ZMULT] = Attributes::makeReal("TMULT", "Multiplier for z dimension.", 1.0);
    itsAttr[Attrib::Distribution::TMULT] = Attributes::makeReal("TMULT", "Multiplier for t dimension.", 1.0);

    itsAttr[Attrib::Distribution::PXMULT] = Attributes::makeReal("PXMULT", "Multiplier for px dimension.", 1.0);
    itsAttr[Attrib::Distribution::PYMULT] = Attributes::makeReal("PYMULT", "Multiplier for py dimension.", 1.0);
    itsAttr[Attrib::Distribution::PZMULT] = Attributes::makeReal("PZMULT", "Multiplier for pz dimension.", 1.0);

    itsAttr[Attrib::Distribution::OFFSETX]
        = Attributes::makeReal("OFFSETX", "Average x offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETY]
        = Attributes::makeReal("OFFSETY", "Average y offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETZ]
        = Attributes::makeReal("OFFSETZ", "Average z offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETT]
        = Attributes::makeReal("OFFSETT", "Average t offset of distribution.", 0.0);

    itsAttr[Attrib::Distribution::OFFSETPX]
        = Attributes::makeReal("OFFSETPX", "Average px offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETPY]
        = Attributes::makeReal("OFFSETPY", "Average py offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETPZ]
        = Attributes::makeReal("OFFSETPZ", "Average pz offset of distribution.", 0.0);
    itsAttr[Attrib::Distribution::OFFSETP]
        = Attributes::makeReal("OFFSETP", "Momentum shift relative to referenc particle.", 0.0);

    // Parameters for defining an initial distribution.
    itsAttr[Attrib::Distribution::SIGMAX] = Attributes::makeReal("SIGMAX", "SIGMAx (m)", 0.0);
    itsAttr[Attrib::Distribution::SIGMAY] = Attributes::makeReal("SIGMAY", "SIGMAy (m)", 0.0);
    itsAttr[Attrib::Distribution::SIGMAR] = Attributes::makeReal("SIGMAR", "SIGMAr (m)", 0.0);
    itsAttr[Attrib::Distribution::SIGMAZ] = Attributes::makeReal("SIGMAZ", "SIGMAz (m)", 0.0);
    itsAttr[Attrib::Distribution::SIGMAT] = Attributes::makeReal("SIGMAT", "SIGMAt (m)", 0.0);
    itsAttr[Attrib::Distribution::TPULSEFWHM]
        = Attributes::makeReal("TPULSEFWHM", "Pulse FWHM for emitted distribution.", 0.0);
    itsAttr[Attrib::Distribution::TRISE]
        = Attributes::makeReal("TRISE", "Rise time for emitted distribution.", 0.0);
    itsAttr[Attrib::Distribution::TFALL]
        = Attributes::makeReal("TFALL", "Fall time for emitted distribution.", 0.0);
    itsAttr[Attrib::Distribution::SIGMAPX] = Attributes::makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[Attrib::Distribution::SIGMAPY] = Attributes::makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[Attrib::Distribution::SIGMAPZ] = Attributes::makeReal("SIGMAPZ", "SIGMApz", 0.0);

    itsAttr[Attrib::Distribution::MX]
        = Attributes::makeReal("MX", "Defines the binomial distribution in x, "
                               "0.0 ... infinity.", 10001.0);
    itsAttr[Attrib::Distribution::MY]
        = Attributes::makeReal("MY", "Defines the binomial distribution in y, "
                               "0.0 ... infinity.", 10001.0);
    itsAttr[Attrib::Distribution::MZ]
        = Attributes::makeReal("MZ", "Defines the binomial distribution in z, "
                               "0.0 ... infinity.", 10001.0);
    itsAttr[Attrib::Distribution::MT]
        = Attributes::makeReal("MT", "Defines the binomial distribution in t, "
                               "0.0 ... infinity.", 10001.0);

    itsAttr[Attrib::Distribution::CUTOFFX] = Attributes::makeReal("CUTOFFX", "Distribution cutoff x "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFY] = Attributes::makeReal("CUTOFFY", "Distribution cutoff r "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFR] = Attributes::makeReal("CUTOFFR", "Distribution cutoff radial "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFLONG]
        = Attributes::makeReal("CUTOFFLONG", "Distribution cutoff z or t direction in "
                               "units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFPX] = Attributes::makeReal("CUTOFFPX", "Distribution cutoff px "
                                                          "dimension in units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFPY] = Attributes::makeReal("CUTOFFPY", "Distribution cutoff py "
                                                          "dimension in units of sigma.", 3.0);
    itsAttr[Attrib::Distribution::CUTOFFPZ] = Attributes::makeReal("CUTOFFPZ", "Distribution cutoff pz "
                                                          "dimension in units of sigma.", 3.0);

    itsAttr[Attrib::Distribution::FTOSCAMPLITUDE]
        = Attributes::makeReal("FTOSCAMPLITUDE", "Amplitude of oscillations superimposed "
                               "on flat top portion of emitted GAUSS "
                               "distribtuion (in percent of flat top "
                               "amplitude)",0.0);
    itsAttr[Attrib::Distribution::FTOSCPERIODS]
        = Attributes::makeReal("FTOSCPERIODS", "Number of oscillations superimposed on "
                               "flat top portion of emitted GAUSS "
                               "distribution", 0.0);

    /*
     * TODO: Find out what these correlations really mean and write
     * good descriptions for them.
     */
    itsAttr[Attrib::Distribution::CORRX]
        = Attributes::makeReal("CORRX", "x/px correlation, (R12 in transport "
                               "notation).", 0.0);
    itsAttr[Attrib::Distribution::CORRY]
        = Attributes::makeReal("CORRY", "y/py correlation, (R34 in transport "
                               "notation).", 0.0);
    itsAttr[Attrib::Distribution::CORRZ]
        = Attributes::makeReal("CORRZ", "z/pz correlation, (R56 in transport "
                               "notation).", 0.0);
    itsAttr[Attrib::Distribution::CORRT]
        = Attributes::makeReal("CORRT", "t/pt correlation, (R56 in transport "
                               "notation).", 0.0);

    itsAttr[Attrib::Distribution::R51]
        = Attributes::makeReal("R51", "x/z correlation, (R51 in transport "
                               "notation).", 0.0);
    itsAttr[Attrib::Distribution::R52]
        = Attributes::makeReal("R52", "xp/z correlation, (R52 in transport "
                               "notation).", 0.0);

    itsAttr[Attrib::Distribution::R61]
        = Attributes::makeReal("R61", "x/pz correlation, (R61 in transport "
                               "notation).", 0.0);
    itsAttr[Attrib::Distribution::R62]
        = Attributes::makeReal("R62", "xp/pz correlation, (R62 in transport "
                               "notation).", 0.0);

    itsAttr[Attrib::Distribution::R]
        = Attributes::makeRealArray("R", "r correlation");

    // Parameters for using laser profile to generate a distribution.
    itsAttr[Attrib::Distribution::LASERPROFFN]
        = Attributes::makeString("LASERPROFFN", "File containing a measured laser image "
                                 "profile (x,y).", "");
    itsAttr[Attrib::Distribution::IMAGENAME]
        = Attributes::makeString("IMAGENAME", "Name of the laser image.", "");
    itsAttr[Attrib::Distribution::INTENSITYCUT]
        = Attributes::makeReal("INTENSITYCUT", "For background subtraction of laser "
                               "image profile, in % of max intensity.",
                               0.0);
    itsAttr[Attrib::Distribution::FLIPX]
        = Attributes::makeBool("FLIPX", "Flip laser profile horizontally", false);
    itsAttr[Attrib::Distribution::FLIPY]
        = Attributes::makeBool("FLIPY", "Flip laser profile vertically", false);
    itsAttr[Attrib::Distribution::ROTATE90]
        = Attributes::makeBool("ROTATE90", "Rotate laser profile 90 degrees counter clockwise", false);
    itsAttr[Attrib::Distribution::ROTATE180]
        = Attributes::makeBool("ROTATE180", "Rotate laser profile 180 degrees counter clockwise", false);
    itsAttr[Attrib::Distribution::ROTATE270]
        = Attributes::makeBool("ROTATE270", "Rotate laser profile 270 degrees counter clockwise", false);

    // Dark current and field emission model parameters.
    itsAttr[Attrib::Distribution::NPDARKCUR]
        = Attributes::makeReal("NPDARKCUR", "Number of dark current particles.", 1000.0);
    itsAttr[Attrib::Distribution::INWARDMARGIN]
        = Attributes::makeReal("INWARDMARGIN", "Inward margin of initialized dark "
                               "current particle positions.", 0.001);
    itsAttr[Attrib::Distribution::EINITHR]
        = Attributes::makeReal("EINITHR", "E field threshold (MV/m), only in position r "
                               "with E(r)>EINITHR that particles will be "
                               "initialized.", 0.0);
    itsAttr[Attrib::Distribution::FNA]
        = Attributes::makeReal("FNA", "Empirical constant A for Fowler-Nordheim "
                               "emission model.", 1.54e-6);
    itsAttr[Attrib::Distribution::FNB]
        = Attributes::makeReal("FNB", "Empirical constant B for Fowler-Nordheim "
                               "emission model.", 6.83e9);
    itsAttr[Attrib::Distribution::FNY]
        = Attributes::makeReal("FNY", "Constant for image charge effect parameter y(E) "
                               "in Fowler-Nordheim emission model.", 3.795e-5);
    itsAttr[Attrib::Distribution::FNVYZERO]
        = Attributes::makeReal("FNVYZERO", "Zero order constant for v(y) function in "
                               "Fowler-Nordheim emission model.", 0.9632);
    itsAttr[Attrib::Distribution::FNVYSECOND]
        = Attributes::makeReal("FNVYSECOND", "Second order constant for v(y) function "
                               "in Fowler-Nordheim emission model.", 1.065);
    itsAttr[Attrib::Distribution::FNPHIW]
        = Attributes::makeReal("FNPHIW", "Work function of gun surface material (eV).",
                               4.65);
    itsAttr[Attrib::Distribution::FNBETA]
        = Attributes::makeReal("FNBETA", "Field enhancement factor for Fowler-Nordheim "
                               "emission.", 50.0);
    itsAttr[Attrib::Distribution::FNFIELDTHR]
        = Attributes::makeReal("FNFIELDTHR", "Field threshold for Fowler-Nordheim "
                               "emission (MV/m).", 30.0);
    itsAttr[Attrib::Distribution::FNMAXEMI]
        = Attributes::makeReal("FNMAXEMI", "Maximum number of electrons emitted from a "
                               "single triangle for Fowler-Nordheim "
                               "emission.", 20.0);
    itsAttr[Attrib::Distribution::SECONDARYFLAG]
        = Attributes::makeReal("SECONDARYFLAG", "Select the secondary model type 0:no "
                               "secondary emission; 1:Furman-Pivi; "
                               "2 or larger: Vaughan's model.", 0.0);
    itsAttr[Attrib::Distribution::NEMISSIONMODE]
        = Attributes::makeBool("NEMISSIONMODE", "Secondary emission mode type true: "
                               "emit n true secondaries; false: emit "
                               "one particle with n times charge.", true);
    itsAttr[Attrib::Distribution::VSEYZERO]
        = Attributes::makeReal("VSEYZERO", "Sey_0 in Vaughan's model.", 0.5);
    itsAttr[Attrib::Distribution::VEZERO]
        = Attributes::makeReal("VEZERO", "Energy related to sey_0 in Vaughan's model "
                               "in eV.", 12.5);
    itsAttr[Attrib::Distribution::VSEYMAX]
        = Attributes::makeReal("VSEYMAX", "Sey max in Vaughan's model.", 2.22);
    itsAttr[Attrib::Distribution::VEMAX]
        = Attributes::makeReal("VEMAX", "Emax in Vaughan's model in eV.", 165.0);
    itsAttr[Attrib::Distribution::VKENERGY]
        = Attributes::makeReal("VKENERGY", "Fitting parameter denotes the roughness of "
                               "surface for impact energy in Vaughan's "
                               "model.", 1.0);
    itsAttr[Attrib::Distribution::VKTHETA]
        = Attributes::makeReal("VKTHETA", "Fitting parameter denotes the roughness of "
                               "surface for impact angle in Vaughan's "
                               "model.", 1.0);
    itsAttr[Attrib::Distribution::VVTHERMAL]
        = Attributes::makeReal("VVTHERMAL", "Thermal velocity of Maxwellian distribution "
                               "of secondaries in Vaughan's model.", 7.268929821 * 1e5);
    itsAttr[Attrib::Distribution::VW]
        = Attributes::makeReal("VW", "VW denote the velocity scalar for parallel plate "
                               "benchmark.", 1.0);
    itsAttr[Attrib::Distribution::SURFMATERIAL]
        = Attributes::makeReal("SURFMATERIAL", "Material type number of the cavity "
                               "surface for Furman-Pivi's model, 0 "
                               "for copper, 1 for stainless steel.", 0.0);

    /*
     *   Feature request Issue #14
     */

    itsAttr[Attrib::Distribution::ID1]
        = Attributes::makeRealArray("ID1", "User defined particle with ID=1");
    itsAttr[Attrib::Distribution::ID2]
        = Attributes::makeRealArray("ID2", "User defined particle with ID=2");


    itsAttr[Attrib::Distribution::SCALABLE]
        = Attributes::makeBool("SCALABLE", "If true then distribution is scalable with "
                               "respect of number of particles and number of cores", false);

    /*
     * Legacy attributes (or ones that need to be implemented.)
     */

    // Parameters for emitting a distribution.
    // itsAttr[Attrib::Legacy::Distribution::DEBIN]
    //     = Attributes::makeReal("DEBIN", "Energy band for a bin in keV that defines "
    //                            "when to combine bins. That is, when the energy "
    //                            "spread of all bins is below this number "
    //                            "combine bins into a single bin.", 1000000.0);
    itsAttr[Attrib::Legacy::Distribution::SBIN]
        = Attributes::makeReal("SBIN", "Number of sample bins to use per energy bin "
                               "when emitting a distribution.", 100.0);
    /*
     * Specific to type GAUSS and GUNGAUSSFLATTOPTH and ASTRAFLATTOPTH. These
     * last two distribution will eventually just become a subset of GAUSS.
     */
    itsAttr[Attrib::Legacy::Distribution::SIGMAPT] = Attributes::makeReal("SIGMAPT", "SIGMApt", 0.0);

    itsAttr[Attrib::Legacy::Distribution::CUTOFF]
        = Attributes::makeReal("CUTOFF", "Longitudinal cutoff for Gaussian in units "
                               "of sigma.", 3.0);


    // Mixed use attributes (used by more than one distribution type).
    itsAttr[Attrib::Legacy::Distribution::T]
        = Attributes::makeReal("T", "Not supported anymore");

    itsAttr[Attrib::Legacy::Distribution::PT]
        = Attributes::makeReal("PT", "Not supported anymore.");


    // Attributes that are not yet implemented.
    // itsAttr[Attrib::Legacy::Distribution::ALPHAX]
    //     = Attributes::makeReal("ALPHAX", "Courant Snyder parameter.", 0.0);
    // itsAttr[Attrib::Legacy::Distribution::ALPHAY]
    //     = Attributes::makeReal("ALPHAY", "Courant Snyder parameter.", 0.0);
    // itsAttr[Attrib::Legacy::Distribution::BETAX]
    //     = Attributes::makeReal("BETAX", "Courant Snyder parameter.", 1.0);
    // itsAttr[Attrib::Legacy::Distribution::BETAY]
    //     = Attributes::makeReal("BETAY", "Courant Snyder parameter.", 1.0);

    // itsAttr[Attrib::Legacy::Distribution::DX]
    //     = Attributes::makeReal("DX", "Dispersion in x (R16 in Transport notation).", 0.0);
    // itsAttr[Attrib::Legacy::Distribution::DDX]
    //     = Attributes::makeReal("DDX", "First derivative of DX.", 0.0);

    // itsAttr[Attrib::Legacy::Distribution::DY]
    //     = Attributes::makeReal("DY", "Dispersion in y (R36 in Transport notation).", 0.0);
    // itsAttr[Attrib::Legacy::Distribution::DDY]
    //     = Attributes::makeReal("DDY", "First derivative of DY.", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

void Distribution::setFieldEmissionParameters() {

    darkCurrentParts_m = static_cast<size_t> (Attributes::getReal(itsAttr[Attrib::Distribution::NPDARKCUR]));
    darkInwardMargin_m = Attributes::getReal(itsAttr[Attrib::Distribution::INWARDMARGIN]);
    eInitThreshold_m   = Attributes::getReal(itsAttr[Attrib::Distribution::EINITHR]);
    workFunction_m     = Attributes::getReal(itsAttr[Attrib::Distribution::FNPHIW]);
    fieldEnhancement_m = Attributes::getReal(itsAttr[Attrib::Distribution::FNBETA]);
    maxFN_m            = static_cast<size_t> (Attributes::getReal(itsAttr[Attrib::Distribution::FNMAXEMI]));
    fieldThrFN_m       = Attributes::getReal(itsAttr[Attrib::Distribution::FNFIELDTHR]);
    paraFNA_m          = Attributes::getReal(itsAttr[Attrib::Distribution::FNA]);
    paraFNB_m          = Attributes::getReal(itsAttr[Attrib::Distribution::FNB]);
    paraFNY_m          = Attributes::getReal(itsAttr[Attrib::Distribution::FNY]);
    paraFNVYZe_m       = Attributes::getReal(itsAttr[Attrib::Distribution::FNVYZERO]);
    paraFNVYSe_m       = Attributes::getReal(itsAttr[Attrib::Distribution::FNVYSECOND]);
    secondaryFlag_m    = Attributes::getReal(itsAttr[Attrib::Distribution::SECONDARYFLAG]);
    ppVw_m             = Attributes::getReal(itsAttr[Attrib::Distribution::VW]);
    vVThermal_m        = Attributes::getReal(itsAttr[Attrib::Distribution::VVTHERMAL]);

}

void Distribution::setDistToEmitted(bool emitted) {
    emitting_m = emitted;
}

void Distribution::setDistType() {
    if (itsAttr[Attrib::Legacy::Distribution::DISTRIBUTION]) {
        throw OpalException("Distribution::setDistType()",
                            "The attribute DISTRIBUTION isn't supported any more, use TYPE instead");
    }

    distT_m = Util::toUpper(Attributes::getString(itsAttr[Attrib::Distribution::TYPE]));
    if (distT_m == "FROMFILE")
        distrTypeT_m = DistrTypeT::FROMFILE;
    else if (distT_m == "GAUSS")
        distrTypeT_m = DistrTypeT::GAUSS;
    else if (distT_m == "BINOMIAL")
        distrTypeT_m = DistrTypeT::BINOMIAL;
    else if (distT_m == "FLATTOP")
        distrTypeT_m = DistrTypeT::FLATTOP;
    else if (distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = DistrTypeT::GUNGAUSSFLATTOPTH;
    else if (distT_m == "ASTRAFLATTOPTH")
        distrTypeT_m = DistrTypeT::ASTRAFLATTOPTH;
    else if (distT_m == "SURFACEEMISSION")
        distrTypeT_m = DistrTypeT::SURFACEEMISSION;
    else if (distT_m == "SURFACERANDCREATE")
        distrTypeT_m = DistrTypeT::SURFACERANDCREATE;
    else if (distT_m == "GAUSSMATCHED")
        distrTypeT_m = DistrTypeT::MATCHEDGAUSS;
    else {
        throw OpalException("Distribution::setDistType()",
                            "The distribution \"" + distT_m + "\" isn't known.\n" +
                            "Known distributions are:\n"
                            "FROMFILE\n"
                            "GAUSS\n"
                            "BINOMIAL\n"
                            "FLATTOP\n"
                            "GUNGAUSSFLATTOPTH\n"
                            "ASTRAFLATTTOPTH\n"
                            "SURFACEEMISSION\n"
                            "SURFACERANDCREATE\n"
                            "GAUSSMATCHED");
    }
}

void Distribution::setEmissionTime(double &maxT, double &minT) {

    if (addedDistributions_m.size() == 0) {

        switch (distrTypeT_m) {

        case DistrTypeT::FLATTOP:
        case DistrTypeT::GAUSS:
        case DistrTypeT::GUNGAUSSFLATTOPTH:
            tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
                * (sigmaTRise_m + sigmaTFall_m);
            break;

        case DistrTypeT::ASTRAFLATTOPTH:
            /*
             * Don't do anything. Emission time is set during the distribution
             * creation. Only this distribution type does it this way. This is
             * a legacy feature.
             */
            break;

        default:
            /*
             * Increase maxT and decrease minT by 0.05% of total time
             * to ensure that no particles fall on the leading edge of
             * the first emission time step or the trailing edge of the last
             * emission time step.
             */
            double deltaT = maxT - minT;
            maxT += deltaT * 0.0005;
            minT -= deltaT * 0.0005;
            tEmission_m = (maxT - minT);
            break;
        }

    } else {
        /*
         * Increase maxT and decrease minT by 0.05% of total time
         * to ensure that no particles fall on the leading edge of
         * the first emission time step or the trailing edge of the last
         * emission time step.
         */
        double deltaT = maxT - minT;
        maxT += deltaT * 0.0005;
        minT -= deltaT * 0.0005;
        tEmission_m = (maxT - minT);
    }
    tBin_m = tEmission_m / numberOfEnergyBins_m;
}

void Distribution::setDistParametersBinomial(double massIneV) {

    /*
     * Set Distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     */
    std::vector<double> cr = Attributes::getRealArray(itsAttr[Attrib::Distribution::R]);

    if (cr.size()>0) {
        throw OpalException("Distribution::setDistParametersBinomial",
                            "Attribute R is not supported for binomial distribution\n"
                            "use CORR[X|Y|Z] and R51, R52, R61, R62 instead");
    }

    correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRX]);
    correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRY]);
    correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRT]);
    correlationMatrix_m(4, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::R51]);
    correlationMatrix_m(4, 1) = Attributes::getReal(itsAttr[Attrib::Distribution::R52]);
    correlationMatrix_m(5, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::R61]);
    correlationMatrix_m(5, 1) = Attributes::getReal(itsAttr[Attrib::Distribution::R62]);

    //CORRZ overrides CORRT. We initially use CORRT for legacy compatibility.
    if (!itsAttr[Attrib::Distribution::CORRZ].defaultUsed())
        correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]);

    sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAX])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAY])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT])));

    // SIGMAZ overrides SIGMAT. We initially use SIGMAT for legacy compatibility.
    if (Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAZ]) != 0.0)
        sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAZ]));

    sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPX])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPY])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPZ])));

    // SIGMAPZ overrides SIGMAPT.
    if (!itsAttr[Attrib::Legacy::Distribution::SIGMAPT].defaultUsed()) {
        WARNMSG("The attribute SIGMAPT may be removed in a future version\n"
                << "use  SIGMAPZ instead" << endl;)
        sigmaP_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Legacy::Distribution::SIGMAPT]));
    }
    // Check what input units we are using for momentum.
    if (inputMoUnits_m == InputMomentumUnitsT::EV) {
        sigmaP_m[0] = converteVToBetaGamma(sigmaP_m[0], massIneV);
        sigmaP_m[1] = converteVToBetaGamma(sigmaP_m[1], massIneV);
        sigmaP_m[2] = converteVToBetaGamma(sigmaP_m[2], massIneV);
    }

    mBinomial_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::MX])),
                           std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::MY])),
                           std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::MT])));

    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFLONG]));

    cutoffP_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPZ]));

    if (mBinomial_m[2] == 0.0
        || Attributes::getReal(itsAttr[Attrib::Distribution::MZ]) != 0.0)
        mBinomial_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::MZ]));

    if (emitting_m) {
        sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));
        mBinomial_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::MT]));
        correlationMatrix_m(5, 4) = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::CORRT]));
    }
}

void Distribution::setDistParametersFlattop(double massIneV) {

    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     */
    sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPX])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPY])),
                        std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPZ])));

    // Check what input units we are using for momentum.
    if (inputMoUnits_m == InputMomentumUnitsT::EV) {
        sigmaP_m[0] = converteVToBetaGamma(sigmaP_m[0], massIneV);
        sigmaP_m[1] = converteVToBetaGamma(sigmaP_m[1], massIneV);
        sigmaP_m[2] = converteVToBetaGamma(sigmaP_m[2], massIneV);
    }

    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFLONG]));

    correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRX]);
    correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRY]);
    correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRT]);

    // CORRZ overrides CORRT.
    if (Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]) != 0.0)
        correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]);


    if (emitting_m) {
        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAY])),
                            0.0);

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));

        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TPULSEFWHM]));

        /*
         * If TRISE and TFALL are defined > 0.0 then these attributes
         * override SIGMAT.
         */
        if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0) {

            double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL]))
                    / timeRatio;

        }

        // For an emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);

    } else {

        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAY])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAZ])));

    }

    if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR])) > 0.0) {
        sigmaR_m[0] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR]));
        sigmaR_m[1] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR]));
    }

    // Set laser profile/
    laserProfileFileName_m = Attributes::getString(itsAttr[Attrib::Distribution::LASERPROFFN]);
    if (!(laserProfileFileName_m == std::string(""))) {
        laserImageName_m = Attributes::getString(itsAttr[Attrib::Distribution::IMAGENAME]);
        laserIntensityCut_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::INTENSITYCUT]));
        short flags = 0;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::FLIPX])) flags |= LaserProfile::FLIPX;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::FLIPY])) flags |= LaserProfile::FLIPY;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE90])) flags |= LaserProfile::ROTATE90;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE180])) flags |= LaserProfile::ROTATE180;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE270])) flags |= LaserProfile::ROTATE270;

        laserProfile_m = new LaserProfile(laserProfileFileName_m,
                                          laserImageName_m,
                                          laserIntensityCut_m,
                                          flags);
    }

    // Legacy for ASTRAFLATTOPTH.
    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        tRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]));

}

void Distribution::setDistParametersGauss(double massIneV) {

    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     * In case of DistrTypeT::MATCHEDGAUSS we only need to set the cutoff parameters
     */


    cutoffP_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPZ]));


    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFLONG]));

    if  (distrTypeT_m != DistrTypeT::MATCHEDGAUSS) {
        sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPX])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPY])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPZ])));

        // SIGMAPZ overrides SIGMAPT. We initially use SIGMAPT for legacy compatibility.
        if (Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPZ]) != 0.0)
            sigmaP_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAPZ]));

        // Check what input units we are using for momentum.
        if (inputMoUnits_m == InputMomentumUnitsT::EV) {
            sigmaP_m[0] = converteVToBetaGamma(sigmaP_m[0], massIneV);
            sigmaP_m[1] = converteVToBetaGamma(sigmaP_m[1], massIneV);
            sigmaP_m[2] = converteVToBetaGamma(sigmaP_m[2], massIneV);
        }

        std::vector<double> cr = Attributes::getRealArray(itsAttr[Attrib::Distribution::R]);

        if (cr.size()>0) {
            if (cr.size() == 15) {
                *gmsg << "* Use r to specify correlations" << endl;
                unsigned int k = 0;
                for (unsigned int i = 0; i < 5; ++ i) {
                    for (unsigned int j = i + 1; j < 6; ++ j, ++ k) {
                        correlationMatrix_m(j, i) = cr.at(k);
                    }
                }

            }
            else {
                throw OpalException("Distribution::SetDistParametersGauss",
                                    "Inconsistent set of correlations specified, check manual");
            }
        }
        else {
            correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRX]);
            correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRY]);
            correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRT]);
            correlationMatrix_m(4, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::R51]);
            correlationMatrix_m(4, 1) = Attributes::getReal(itsAttr[Attrib::Distribution::R52]);
            correlationMatrix_m(5, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::R61]);
            correlationMatrix_m(5, 1) = Attributes::getReal(itsAttr[Attrib::Distribution::R62]);

            // CORRZ overrides CORRT.
            if (Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]) != 0.0)
                correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]);
        }
    }

    if (emitting_m) {

        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAY])),
                            0.0);

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));

        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TPULSEFWHM]));

        /*
         * If TRISE and TFALL are defined then these attributes
         * override SIGMAT.
         */
        if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0) {

            double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL]))
                    / timeRatio;

        }

        // For and emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);

    } else {
        if  (distrTypeT_m != DistrTypeT::MATCHEDGAUSS) {
            sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAX])),
                                std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAY])),
                                std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT])));

            // SIGMAZ overrides SIGMAT. We initially use SIGMAT for legacy compatibility.
            if (Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAZ]) != 0.0)
                sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAZ]));

        }
    }

    if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR])) > 0.0) {
        sigmaR_m[0] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR]));
        sigmaR_m[1] = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR]));
        cutoffR_m[0] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
        cutoffR_m[1] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
    }
}

void Distribution::setupEmissionModel(PartBunchBase<double, 3> *beam) {

    std::string model = Util::toUpper(Attributes::getString(itsAttr[Attrib::Distribution::EMISSIONMODEL]));
    if (model == "ASTRA")
        emissionModel_m = EmissionModelT::ASTRA;
    else if (model == "NONEQUIL")
        emissionModel_m = EmissionModelT::NONEQUIL;
    else
        emissionModel_m = EmissionModelT::NONE;

    /*
     * The ASTRAFLATTOPTH  of GUNGAUSSFLATTOPTH distributions always uses the
     * ASTRA emission model.
     */
    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH
        || distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH)
        emissionModel_m = EmissionModelT::ASTRA;

    switch (emissionModel_m) {

    case EmissionModelT::ASTRA:
        setupEmissionModelAstra(beam);
        break;
    case EmissionModelT::NONEQUIL:
        setupEmissionModelNonEquil();
        break;
    default:
        break;
    }

}

void Distribution::setupEmissionModelAstra(PartBunchBase<double, 3> *beam) {

    double wThermal = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]));
    pTotThermal_m = converteVToBetaGamma(wThermal, beam->getM());
    pmean_m = Vector_t(0.0, 0.0, 0.5 * pTotThermal_m);
}

void Distribution::setupEmissionModelNone(PartBunchBase<double, 3> *beam) {

    double wThermal = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]));
    pTotThermal_m = converteVToBetaGamma(wThermal, beam->getM());
    double avgPz = std::accumulate(pzDist_m.begin(), pzDist_m.end(), 0.0);
    size_t numParticles = pzDist_m.size();
    reduce(avgPz, avgPz, OpAddAssign());
    reduce(numParticles, numParticles, OpAddAssign());
    avgPz /= numParticles;

    pmean_m = Vector_t(0.0, 0.0, pTotThermal_m + avgPz);
}

void Distribution::setupEmissionModelNonEquil() {

    cathodeWorkFunc_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::W]));
    laserEnergy_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::ELASER]));
    cathodeFermiEnergy_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::FE]));
    cathodeTemp_m = Physics::kB * std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::CATHTEMP]));

    /*
     * Upper limit on energy distribution theoretically goes to infinity.
     * Practically we limit to a probability of 1 part in a billion.
     */
    emitEnergyUpperLimit_m = cathodeFermiEnergy_m
        + cathodeTemp_m * log(1.0e9 - 1.0);

    // TODO: get better estimate of pmean
    pmean_m = Vector_t(0, 0, sqrt(pow(0.5 * emitEnergyUpperLimit_m / (Physics::m_e * 1e9) + 1.0, 2) - 1.0));
}

void Distribution::setupEnergyBins(double maxTOrZ, double minTOrZ) {

    energyBinHist_m = gsl_histogram_alloc(numberOfSampleBins_m * numberOfEnergyBins_m);

    if (emitting_m)
        gsl_histogram_set_ranges_uniform(energyBinHist_m, -tEmission_m, 0.0);
    else
        gsl_histogram_set_ranges_uniform(energyBinHist_m, minTOrZ, maxTOrZ);

}

void Distribution::setupParticleBins(double massIneV, PartBunchBase<double, 3> *beam) {

    numberOfEnergyBins_m
        = static_cast<int>(std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::NBIN])));

    if (numberOfEnergyBins_m > 0) {
        delete energyBins_m;

        int sampleBins = static_cast<int>(std::abs(Attributes::getReal(itsAttr[Attrib::Legacy::Distribution::SBIN])));
        energyBins_m = new PartBins(numberOfEnergyBins_m, sampleBins);

        if (!itsAttr[Attrib::Legacy::Distribution::PT].defaultUsed())
            throw OpalException("Distribution::setupParticleBins",
                                "PT is obsolete. The moments of the beam is defined with OFFSETPZ");

        // we get gamma from PC of the beam
        const double pz    = beam->getP()/beam->getM();
        double gamma = sqrt(pow(pz, 2.0) + 1.0);
        energyBins_m->setGamma(gamma);

    } else {
        energyBins_m = NULL;
    }
}

void Distribution::shiftBeam(double &maxTOrZ, double &minTOrZ) {

    if (emitting_m) {

        if (addedDistributions_m.size() == 0) {

            if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH) {
                for (double& tOrZ : tOrZDist_m)
                    tOrZ -= tEmission_m / 2.0;

                minTOrZ -= tEmission_m / 2.0;
                maxTOrZ -= tEmission_m / 2.0;
            } else if (distrTypeT_m == DistrTypeT::GAUSS
                       || distrTypeT_m == DistrTypeT::FLATTOP
                       || distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH) {
                for (double& tOrZ : tOrZDist_m)
                    tOrZ -= tEmission_m;

                minTOrZ -= tEmission_m;
                maxTOrZ -= tEmission_m;
            } else {
                for (double& tOrZ : tOrZDist_m)
                    tOrZ -= maxTOrZ;

                minTOrZ -= maxTOrZ;
                maxTOrZ -= maxTOrZ;
            }

        } else {
            for (double& tOrZ : tOrZDist_m)
                tOrZ -= maxTOrZ;

            minTOrZ -= maxTOrZ;
            maxTOrZ -= maxTOrZ;
        }

    } else {
        double avgZ[] = {0.0, 1.0 * tOrZDist_m.size()};
        for (double tOrZ : tOrZDist_m)
            avgZ[0] += tOrZ;

        reduce(avgZ, avgZ + 2, avgZ, OpAddAssign());
        avgZ[0] /= avgZ[1];

        for (double& tOrZ : tOrZDist_m)
            tOrZ -= avgZ[0];
    }
}

double Distribution::getEmissionTimeShift() const {
    if (emitting_m)
        return Attributes::getReal(itsAttr[Attrib::Distribution::OFFSETT]);

    return 0.0;
}

void Distribution::shiftDistCoordinates(double massIneV) {

    size_t startIdx = 0;
    for (unsigned int i = 0; i <= addedDistributions_m.size(); ++ i) {
        Distribution *currDist = this;
        if (i > 0)
            currDist = addedDistributions_m[i - 1];
        double deltaX = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETX]);
        double deltaY = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETY]);

        /*
         * OFFSETZ overrides T if it is nonzero. We initially use T
         * for legacy compatiblity. OFFSETT always overrides T, even
         * when zero, for an emitted beam.
         */
        if (currDist->itsAttr[Attrib::Legacy::Distribution::T]) {
            throw OpalException("Distribution::shiftDistCoordinates",
                                "Attribute T isn't supported anymore; use OFFSETZ instead");
        }

        double deltaTOrZ = 0.0;
        if (!emitting_m)
            if (Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETZ]) != 0.0)
                deltaTOrZ = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETZ]);

        double deltaPx = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETPX]);
        double deltaPy = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETPY]);
        double deltaPz = Attributes::getReal(currDist->itsAttr[Attrib::Distribution::OFFSETPZ]);

        if (Attributes::getReal(currDist->itsAttr[Attrib::Legacy::Distribution::PT])!=0.0)
            WARNMSG("PT & PZ are obsolete and will be ignored. The moments of the beam is defined with PC" << endl);

        // Check input momentum units.
        if (inputMoUnits_m == InputMomentumUnitsT::EV) {
            deltaPx = converteVToBetaGamma(deltaPx, massIneV);
            deltaPy = converteVToBetaGamma(deltaPy, massIneV);
            deltaPz = converteVToBetaGamma(deltaPz, massIneV);
        }

        size_t endIdx = startIdx + particlesPerDist_m[i];
        for (size_t particleIndex = startIdx; particleIndex < endIdx; ++ particleIndex) {
            xDist_m.at(particleIndex) += deltaX;
            pxDist_m.at(particleIndex) += deltaPx;
            yDist_m.at(particleIndex) += deltaY;
            pyDist_m.at(particleIndex) += deltaPy;
            tOrZDist_m.at(particleIndex) += deltaTOrZ;
            pzDist_m.at(particleIndex) += deltaPz;
        }

        startIdx = endIdx;
    }
}

void Distribution::writeOutFileHeader() {

    if (Attributes::getBool(itsAttr[Attrib::Distribution::WRITETOFILE]) == false)
        return;

    unsigned int totalNum = tOrZDist_m.size();
    reduce(totalNum, totalNum, OpAddAssign());
    if (Ippl::myNode() != 0)
        return;

    std::string fname = "data/" + OpalData::getInstance()->getInputBasename() + "_" + getOpalName() + ".dat";
    *gmsg << "\n"
          << std::left << std::setw(84) << std::setfill('*') << "* " << "\n"
          << "* Write initial distribution to file \"" << fname << "\"\n"
          << std::left << std::setw(84) << std::setfill('*') << "* "
          << std::setfill(' ') << endl;

    std::ofstream outputFile(fname);
    if (outputFile.bad()) {
        *gmsg << "Unable to open output file \"" << fname << "\"" << endl;
    } else {
        outputFile.setf(std::ios::left);
        outputFile << "# ";
        if (emitting_m) {
            outputFile.width(17);
            outputFile << "x [m]";
            outputFile.width(17);
            outputFile << "px [betax gamma]";
            outputFile.width(17);
            outputFile << "y [m]";
            outputFile.width(17);
            outputFile << "py [betay gamma]";
            outputFile.width(17);
            outputFile << "t [s]";
            outputFile.width(17);
            outputFile << "pz [betaz gamma]" ;
            outputFile.width(17);
            outputFile << "Bin Number" << std::endl;
        } else {
            outputFile.width(17);
            outputFile << "x [m]";
            outputFile.width(17);
            outputFile << "px [betax gamma]";
            outputFile.width(17);
            outputFile << "y [m]";
            outputFile.width(17);
            outputFile << "py [betay gamma]";
            outputFile.width(17);
            outputFile << "z [m]";
            outputFile.width(17);
            outputFile << "pz [betaz gamma]";
            if (numberOfEnergyBins_m > 0) {
                outputFile.width(17);
                outputFile << "Bin Number";
            }
            outputFile << std::endl;

            outputFile << "# " << totalNum << std::endl;
        }
    }
    outputFile.close();
}

void Distribution::writeOutFileEmission() {

    if (!Attributes::getBool(itsAttr[Attrib::Distribution::WRITETOFILE])) {
        xWrite_m.clear();
        pxWrite_m.clear();
        yWrite_m.clear();
        pyWrite_m.clear();
        tOrZWrite_m.clear();
        pzWrite_m.clear();
        binWrite_m.clear();

        return;
    }

    // Gather particles to be written onto node 0.
    std::vector<char> msgbuf;
    const size_t bitsPerParticle = (6 * sizeof(double) + sizeof(size_t));
    size_t totalSendBits = xWrite_m.size() * bitsPerParticle;

    std::vector<long> numberOfBits(Ippl::getNodes(), 0);
    numberOfBits[Ippl::myNode()] = totalSendBits;

    if (Ippl::myNode() == 0) {
        MPI_Reduce(MPI_IN_PLACE, &(numberOfBits[0]), Ippl::getNodes(), MPI_UNSIGNED_LONG, MPI_SUM, 0, Ippl::getComm());
    } else {
        MPI_Reduce(&(numberOfBits[0]), NULL, Ippl::getNodes(), MPI_UNSIGNED_LONG, MPI_SUM, 0, Ippl::getComm());
    }

    Ippl::Comm->barrier();
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG2, IPPL_APP_CYCLE);
    if (Ippl::myNode() > 0) {
        if (totalSendBits > 0) {
            msgbuf.reserve(totalSendBits);
            const char *buffer;
            for (size_t idx = 0; idx < xWrite_m.size(); ++ idx) {
                buffer = reinterpret_cast<const char*>(&(xWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(pxWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(yWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(pyWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(tOrZWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(pzWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(double));
                buffer = reinterpret_cast<const char*>(&(binWrite_m[idx]));
                msgbuf.insert(msgbuf.end(), buffer, buffer + sizeof(size_t));
            }

            Ippl::Comm->raw_send(&(msgbuf[0]), totalSendBits, 0, tag);
        }
    } else {
        std::string fname = "data/" + OpalData::getInstance()->getInputBasename() + "_" + getOpalName() + ".dat";
        std::ofstream outputFile(fname, std::fstream::app);
        if (outputFile.bad()) {
            ERRORMSG(level1 << "Unable to write to file \"" << fname << "\"" << endl);
            for (int node = 1; node < Ippl::getNodes(); ++ node) {
                if (numberOfBits[node] == 0) continue;
                char *recvbuf = new char[numberOfBits[node]];
                Ippl::Comm->raw_receive(recvbuf, numberOfBits[node], node, tag);
                delete[] recvbuf;
            }
        } else {

            outputFile.precision(9);
            outputFile.setf(std::ios::scientific);
            outputFile.setf(std::ios::right);

            for (size_t partIndex = 0; partIndex < xWrite_m.size(); partIndex++) {

                outputFile << std::setw(17) << xWrite_m.at(partIndex)
                           << std::setw(17) << pxWrite_m.at(partIndex)
                           << std::setw(17) << yWrite_m.at(partIndex)
                           << std::setw(17) << pyWrite_m.at(partIndex)
                           << std::setw(17) << tOrZWrite_m.at(partIndex)
                           << std::setw(17) << pzWrite_m.at(partIndex)
                           << std::setw(17) << binWrite_m.at(partIndex) << std::endl;
            }

            int numSenders = 0;
            for (int i = 1; i < Ippl::getNodes(); ++ i) {
                if (numberOfBits[i] > 0) numSenders ++;
            }

            for (int i = 0; i < numSenders; ++ i) {
                int node = Communicate::COMM_ANY_NODE;
                char *recvbuf;
                const int bufsize = Ippl::Comm->raw_probe_receive(recvbuf, node, tag);

                int j = 0;
                while (j < bufsize) {
                    const double *dbuffer = reinterpret_cast<const double*>(recvbuf + j);
                    const size_t *sbuffer = reinterpret_cast<const size_t*>(recvbuf + j + 6 * sizeof(double));
                    outputFile << std::setw(17) << dbuffer[0]
                               << std::setw(17) << dbuffer[1]
                               << std::setw(17) << dbuffer[2]
                               << std::setw(17) << dbuffer[3]
                               << std::setw(17) << dbuffer[4]
                               << std::setw(17) << dbuffer[5]
                               << std::setw(17) << sbuffer[0]
                               << std::endl;
                    j += bitsPerParticle;
                }

                delete[] recvbuf;

            }
        }
        outputFile.close();

    }

    // Clear write vectors.
    xWrite_m.clear();
    pxWrite_m.clear();
    yWrite_m.clear();
    pyWrite_m.clear();
    tOrZWrite_m.clear();
    pzWrite_m.clear();
    binWrite_m.clear();

}

void Distribution::writeOutFileInjection() {

    if (Attributes::getBool(itsAttr[Attrib::Distribution::WRITETOFILE]) == false)
        return;

    std::string fname = "data/" + OpalData::getInstance()->getInputBasename() + "_" + getOpalName() + ".dat";
    // Nodes take turn writing particles to file.
    for (int nodeIndex = 0; nodeIndex < Ippl::getNodes(); nodeIndex++) {

        // Write to file if its our turn.
        size_t numberOfParticles = 0;
        if (Ippl::myNode() == nodeIndex) {
            std::ofstream outputFile(fname, std::fstream::app);
            if (outputFile.bad()) {
                *gmsg << "Node " << Ippl::myNode() << " unable to write"
                      << "to file \"" << fname << "\"" << endl;
            } else {

                outputFile.precision(9);
                outputFile.setf(std::ios::scientific);
                outputFile.setf(std::ios::right);

                numberOfParticles = tOrZDist_m.size();
                for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

                    outputFile.width(17);
                    outputFile << xDist_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pxDist_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << yDist_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pyDist_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << tOrZDist_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pzDist_m.at(partIndex);
                    if (numberOfEnergyBins_m > 0) {
                        size_t binNumber = findEBin(tOrZDist_m.at(partIndex));
                        outputFile.width(17);
                        outputFile << binNumber;
                    }
                    outputFile << std::endl;

                }
            }
            outputFile.close();
        }

        // Wait for writing node before moving on.
        reduce(numberOfParticles, numberOfParticles, OpAddAssign());
    }
}

double Distribution::MDependentBehavior::get(double rand) {
    return 2.0 * sqrt(1.0 - pow(rand, ami_m));
}

double Distribution::GaussianLikeBehavior::get(double rand) {
    return sqrt(-2.0 * log(rand));
}

void Distribution::adjustPhaseSpace(double massIneV) {
    if (emitting_m || distrTypeT_m == DistrTypeT::FROMFILE || OpalData::getInstance()->isInOPALCyclMode())
        return;

    double deltaPx = Attributes::getReal(itsAttr[Attrib::Distribution::OFFSETPX]);
    double deltaPy = Attributes::getReal(itsAttr[Attrib::Distribution::OFFSETPY]);
    // Check input momentum units.
    if (inputMoUnits_m == InputMomentumUnitsT::EV) {
        deltaPx = converteVToBetaGamma(deltaPx, massIneV);
        deltaPy = converteVToBetaGamma(deltaPy, massIneV);
    }

    double avrg[6];
    avrg[0] = std::accumulate(xDist_m.begin(), xDist_m.end(), 0.0) / totalNumberParticles_m;
    avrg[1] = std::accumulate(pxDist_m.begin(), pxDist_m.end(), 0.0) / totalNumberParticles_m;
    avrg[2] = std::accumulate(yDist_m.begin(), yDist_m.end(), 0.0) / totalNumberParticles_m;
    avrg[3] = std::accumulate(pyDist_m.begin(), pyDist_m.end(), 0.0) / totalNumberParticles_m;
    avrg[4] = std::accumulate(tOrZDist_m.begin(), tOrZDist_m.end(), 0.0) / totalNumberParticles_m;
    avrg[5] = 0.0;
    for (unsigned int i = 0; i < pzDist_m.size(); ++ i) {
        // taylor series of sqrt(px*px + py*py + pz*pz) = pz * sqrt(1 + eps*eps) where eps << 1
        avrg[5] += (pzDist_m[i] +
                    (std::pow(pxDist_m[i] + deltaPx, 2) +
                     std::pow(pyDist_m[i] + deltaPy, 2)) / (2 * pzDist_m[i]));
    }
    allreduce(&(avrg[0]), 6, std::plus<double>());
    avrg[5] /= totalNumberParticles_m;

    // solve
    // \sum_{i = 0}^{N-1} \sqrt{(pz_i + \eps)^2 + px_i^2 + py_i^2} = N p
    double eps = avrgpz_m - avrg[5];
    for (unsigned int i = 0; i < pzDist_m.size(); ++ i) {
        xDist_m[i] -= avrg[0];
        pxDist_m[i] -= avrg[1];
        yDist_m[i] -= avrg[2];
        pyDist_m[i] -= avrg[3];
        tOrZDist_m[i] -= avrg[4];
        pzDist_m[i] += eps;
    }
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: