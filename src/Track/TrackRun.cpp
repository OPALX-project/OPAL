// -----------------------------------------------------------------------
// /*$RCSfile*/: TrackRun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackRun
//   The class for the OPAL RUN command.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackRun.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include "Algorithms/Tracker.h"
#include "Algorithms/ThinTracker.h"
#include "Algorithms/ThickTracker.h"

#include "Algorithms/bet/EnvelopeBunch.h"

#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelSliceTracker.h"
#include "Algorithms/ParallelCyclotronTracker.h"
#include "Algorithms/StatisticalErrors.h"
#include "Algorithms/NilTracker.h"

#include "Attributes/Attributes.h"
#include "Beamlines/TBeamline.h"

#include "BasicActions/Option.h"

#include "Distribution/Distribution.h"
#include "Elements/OpalBeamBeam3D.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include "Utilities/Round.h"
#include "Utilities/Util.h"
#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/FieldSolver.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"
#include "Structure/H5PartWrapperForPC.h"
#include "Structure/H5PartWrapperForPS.h"

#include "OPALconfig.h"
#include "changes.h"

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <iomanip>

extern Inform *gmsg;

using namespace Physics;

// ------------------------------------------------------------------------

namespace {

    // The attributes of class TrackRun.
    enum {
        METHOD,       // Tracking method to use.
        TURNS,        // The number of turns to be tracked.
        MBMODE,       // The working way for multi-bunch mode for OPAL-cycl: "FORCE" or "AUTO"
        PARAMB,       // The control parameter for "AUTO" mode of multi-bunch,
        MB_ETA,       // The scale parameter for binning in multi-bunch mode
        MB_BINNING,   // The binning type in multi-bunch mode
        BEAM,         // The beam to track
        FIELDSOLVER,  // The field solver attached
        BOUNDARYGEOMETRY, // The boundary geometry
        DISTRIBUTION, // The particle distribution
        MULTIPACTING, // MULTIPACTING flag
        OBJECTIVES,   // for which columns of an sdds file statistical errors should be considered
        // THE INTEGRATION TIMESTEP IN SEC
        SIZE
    };
}

const std::string TrackRun::defaultDistribution("DISTRIBUTION");

TrackRun::TrackRun():
    Action(SIZE, "RUN",
           "The \"RUN\" sub-command tracks the defined particles through "
           "the given lattice."),
    itsTracker(NULL),
    dist(NULL),
    fs(NULL),
    ds(NULL),
    phaseSpaceSink_m(NULL) {
    itsAttr[METHOD] = Attributes::makeString
                      ("METHOD", "Name of tracking algorithm to use:\n"
                       "\t\t\t\"THIN\" (default) or \"THICK, OPAL-T,OPAL-T3D, OPAL-CYCL\".", "THIN");
    itsAttr[TURNS] = Attributes::makeReal
                     ("TURNS", "Number of turns to be tracked; Number of neighboring bunches to be tracked in cyclotron", 1.0);

    itsAttr[MBMODE] = Attributes::makeString
                      ("MBMODE", "The working way for multi-bunch mode for OPAL-cycl: FORCE or AUTO ", "FORCE");

    itsAttr[PARAMB] = Attributes::makeReal
                      ("PARAMB", "Control parameter to define when to start multi-bunch mode, only available in \"AUTO\" mode ", 5.0);

    itsAttr[MB_ETA] = Attributes::makeReal("MB_ETA",
                                           "The scale parameter for binning in multi-bunch mode",
                                           0.01);
    
    itsAttr[MB_BINNING] = Attributes::makeString
                          ("MB_BINNING", "Type of energy binning in multi-bunch mode: GAMMA or BUNCH", "GAMMA");

    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "Name of beam ", "BEAM");
    itsAttr[FIELDSOLVER] = Attributes::makeString
                           ("FIELDSOLVER", "Field solver to be used ", "FIELDSOLVER");
    itsAttr[BOUNDARYGEOMETRY] = Attributes::makeString
                           ("BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default)", "NONE");
    itsAttr[DISTRIBUTION] = Attributes::makeStringArray
                             ("DISTRIBUTION", "List of particle distributions to be used ");
    itsAttr[MULTIPACTING] = Attributes::makeBool
                            ("MULTIPACTING", "Multipacting flag, default: false. Set true to initialize primary particles according to BoundaryGeometry", false);
    itsAttr[OBJECTIVES] = Attributes::makeStringArray
                          ("OBJECTIVES", "List of SDDS columns that should be considered when evaluating statistical errors");

    registerOwnership(AttributeHandler::SUB_COMMAND);

    opal = OpalData::getInstance();
}


TrackRun::TrackRun(const std::string &name, TrackRun *parent):
    Action(name, parent),
    itsTracker(NULL),
    dist(NULL),
    fs(NULL),
    ds(NULL),
    phaseSpaceSink_m(NULL) {
    opal = OpalData::getInstance();
}


TrackRun::~TrackRun()
{
    delete phaseSpaceSink_m;
}


TrackRun *TrackRun::clone(const std::string &name) {
    return new TrackRun(name, this);
}


void TrackRun::execute() {
    const int currentVersion = ((OPAL_VERSION_MAJOR * 100) + OPAL_VERSION_MINOR) * 100;
    if (Options::version < currentVersion) {
        unsigned int fileVersion = Options::version / 100;
        bool newerChanges = false;
        for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++ it) {
            if (it->first > fileVersion) {
                newerChanges = true;
                break;
            }
        }
        if (newerChanges) {
            Inform errorMsg("Error");
            errorMsg << "\n******************** V E R S I O N   M I S M A T C H ***********************\n" << endl;
            for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++ it) {
                if (it->first > fileVersion) {
                    errorMsg << it->second << endl;
                }
            }
            errorMsg << "\n"
                     << "* Make sure you do understand these changes and adjust your input file \n"
                     << "* accordingly. Then add\n"
                     << "* OPTION, VERSION = " << currentVersion << ";\n"
                     << "* to your input file. " << endl;
            errorMsg << "\n****************************************************************************\n" << endl;
            throw OpalException("TrackRun::execute", "Version mismatch");
        }
    }

    // Get algorithm to use.
    std::string method = Util::toUpper(Attributes::getString(itsAttr[METHOD]));
    if(method == "THIN") {
        *gmsg << "  Method == \"THIN\"" << endl;
        itsTracker = new ThinTracker(*Track::block->use->fetchLine(),
                                     Track::block->bunch, Track::block->reference,
                                     false, false);
    } else if(method == "THICK") {
	setupThickTracker();
    } else if(method == "PARALLEL-SLICE" || method == "OPAL-E") {
        setupSliceTracker();
    } else if(method == "PARALLEL-T" || method == "OPAL-T") {
        setupTTracker();
    } else if(method == "CYCLOTRON-T" || method == "OPAL-CYCL") {
        setupCyclotronTracker();
    } else if(method.substr(0,18) == "STATISTICAL-ERRORS") {
        setupStatisticalErrors(method);
    } else {
        throw OpalException("TrackRun::execute()",
                            "Method name \"" + method + "\" unknown.");
    }

    if(method == "THIN" || method == "THICK") {
        int turns = int(Round(Attributes::getReal(itsAttr[TURNS])));

        // Track for the all but last turn.
        for(int turn = 1; turn < turns; ++turn) {
            itsTracker->execute();
        }

        // Track the last turn.
        itsTracker->execute();

    } else {
        itsTracker->execute();

        opal->setRestartRun(false);
    }

    opal->bunchIsAllocated();
    if(method == "PARALLEL-SLICE")
        opal->slbunchIsAllocated();

    delete itsTracker;
}

void TrackRun::setupSliceTracker() {
    OpalData::getInstance()->setInOPALEnvMode();
    bool isFollowupTrack = opal->hasSLBunchAllocated();

    if(!opal->hasSLBunchAllocated()) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-SLICE, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else if(isFollowupTrack) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-SLICE, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam   *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    Track::block->bunch->setBeamFrequency(beam->getFrequency() * 1e6);

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack) {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    const size_t numberOfDistributions = distr_str.size();
    if (numberOfDistributions == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distr_str.at(0));
        // dist->setNumberOfDistributions(numberOfDistributions);

        // if(numberOfDistributions > 1) {
        //     *gmsg << "Found more than one distribution: ";
        //     for(size_t i = 1; i < numberOfDistributions; ++ i) {
        //         Distribution *d = Distribution::find(distr_str.at(i));

        //         d->setNumberOfDistributions(numberOfDistributions);
        //         distrs_m.push_back(d);

        //         *gmsg << " " << distr_str.at(i);
        //     }
        //     *gmsg << endl;
        // }
    }

    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
    fs->initCartesianFields();

    double charge = 0.0;

    if(!opal->hasSLBunchAllocated()) {
        if(!opal->inRestartRun()) {

            dist->createOpalE(beam, distrs_m, Track::block->slbunch, 0.0, 0.0);
            opal->setGlobalPhaseShift(0.5 * dist->getTEmission());

        } else {
            /***
                reload slice distribution
            */

            dist->doRestartOpalE(Track::block->slbunch, beam->getNumberOfParticles(), opal->getRestartStep(), phaseSpaceSink_m);
        }
    } else {
        charge = 1.0;
    }

    Track::block->slbunch->setdT(Track::block->dT.front());
    // set the total charge
    charge = beam->getCharge() * beam->getCurrent() / (beam->getFrequency() * 1e6);
    Track::block->slbunch->setCharge(charge);
    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->slbunch->setCouplingConstant(coefE);
    //Track::block->slbunch->calcBeamParameters();


    if(!opal->inRestartRun()) {
        if(!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, -1));
    }

    ds = opal->getDataSink();

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = Track::block->zstart;
        auto &zstop = Track::block->zstop;
        auto &timeStep = Track::block->localTimeSteps;
        auto &dT = Track::block->dT;

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
        }

        zstop.erase(zstop.begin(), zstop.begin() + i);
        timeStep.erase(timeStep.begin(), timeStep.begin() + i);
        dT.erase(dT.begin(), dT.begin() + i);

        Track::block->bunch->setdT(dT.front());
    } else {
        Track::block->zstart = 0.0;
    }

    if(!opal->hasBunchAllocated())
        *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *Track::block->slbunch  << endl;
    *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq
          << " Inputfile is " << opal->getInputFn() << endl;


    // findPhasesForMaxEnergy();

    itsTracker = new ParallelSliceTracker(*Track::block->use->fetchLine(),
                                          dynamic_cast<EnvelopeBunch &>(*Track::block->slbunch),
                                          *ds,
                                          Track::block->reference,
                                          false,
                                          false,
                                          Track::block->localTimeSteps,
                                          Track::block->zstart,
                                          Track::block->zstop,
                                          Track::block->dT);
}



void TrackRun::setupThickTracker()
{
    bool isFollowupTrack = opal->hasBunchAllocated();

    if(isFollowupTrack) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == THICK, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
        Track::block->bunch->setLocalTrackStep(0);
    } else {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == THICK, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            BoundaryGeometry *bg = BoundaryGeometry::find(Attributes::getString(itsAttr[BOUNDARYGEOMETRY]))->
                clone(getOpalName() + std::string("_geometry"));
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }

    setupFieldsolver();

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    double charge = setDistributionParallelT(beam);

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack && !opal->inRestartRun())
        Track::block->bunch->setT(Track::block->t0_m);
    if (Track::block->bunch->getIfBeamEmitting()) {
      Track::block->bunch->setChargeZeroPart(charge);
    } else {
      Track::block->bunch->setCharge(charge);
    }
    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);


    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    if(!opal->inRestartRun()) {
        if(!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, -1));//opal->getRestartStep()));
    }

    ds = opal->getDataSink();

    if(!opal->hasBunchAllocated())
      *gmsg << *dist << endl;

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = /*Track::block->bunch->get_sPos() +*/ Track::block->zstart;
        auto &zstop = Track::block->zstop;
        auto &timeStep = Track::block->localTimeSteps;
        auto &dT = Track::block->dT;

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
        }

        zstop.erase(zstop.begin(), zstop.begin() + i);
        timeStep.erase(timeStep.begin(), timeStep.begin() + i);
        dT.erase(dT.begin(), dT.begin() + i);

        Track::block->bunch->setdT(dT.front());
    } else {
        Track::block->zstart = 0.0;
    }

    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    *gmsg << level2
          << "Phase space dump frequency " << Options::psDumpFreq << " and "
          << "statistics dump frequency " << Options::statDumpFreq << " w.r.t. the time step." << endl;

    itsTracker = new ThickTracker(*Track::block->use->fetchLine(),
                                  Track::block->bunch, *beam, *ds, Track::block->reference,
                                  false, false, Track::block->localTimeSteps,
                                  Track::block->zstart, Track::block->zstop, Track::block->dT,
                                  Track::block->truncOrder);
}


void TrackRun::setupTTracker(){
    OpalData::getInstance()->setInOPALTMode();
    bool isFollowupTrack = opal->hasBunchAllocated();

    if(isFollowupTrack) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
        Track::block->bunch->setLocalTrackStep(0);
    } else {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-T, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    Track::block->bunch->setBeamFrequency(beam->getFrequency() * 1e6);

    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            BoundaryGeometry *bg = BoundaryGeometry::find(Attributes::getString(itsAttr[BOUNDARYGEOMETRY]))->
                clone(getOpalName() + std::string("_geometry"));
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }

    setupFieldsolver();

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    double charge = setDistributionParallelT(beam);

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack && !opal->inRestartRun())
        Track::block->bunch->setT(Track::block->t0_m);

    bool mpacflg = Attributes::getBool(itsAttr[MULTIPACTING]);
    if(!mpacflg) {
        if (Track::block->bunch->getIfBeamEmitting()) {
            Track::block->bunch->setChargeZeroPart(charge);
        } else {
            Track::block->bunch->setCharge(charge);
        }
        // set coupling constant
        double coefE = 1.0 / (4 * pi * epsilon_0);
        Track::block->bunch->setCouplingConstant(coefE);


        // statistical data are calculated (rms, eps etc.)
        Track::block->bunch->calcBeamParameters();
    } else {
        Track::block->bunch->setChargeZeroPart(charge);// just set bunch->qi_m=charge, don't set bunch->Q[] as we have not initialized any particle yet.
        Track::block->bunch->calcBeamParametersInitial();// we have not initialized any particle yet.
    }

    if(!opal->inRestartRun()) {
        if(!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, -1));//opal->getRestartStep()));
    }

    ds = opal->getDataSink();

    if(!opal->hasBunchAllocated()) {
        if(!mpacflg) {
            *gmsg << std::scientific;
            *gmsg << *dist << endl;
        } else {
            *gmsg << "* Multipacting flag is true. The particle distribution in the run command will be ignored " << endl;
        }
    }

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = Track::block->zstart;
        auto &zstop = Track::block->zstop;
        auto &timeStep = Track::block->localTimeSteps;
        auto &dT = Track::block->dT;

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
        }

        zstop.erase(zstop.begin(), zstop.begin() + i);
        timeStep.erase(timeStep.begin(), timeStep.begin() + i);
        dT.erase(dT.begin(), dT.begin() + i);

        Track::block->bunch->setdT(dT.front());
    } else {
        Track::block->zstart = 0.0;
    }

    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;

    // findPhasesForMaxEnergy();

    *gmsg << level2
          << "Phase space dump frequency " << Options::psDumpFreq << " and "
          << "statistics dump frequency " << Options::statDumpFreq << " w.r.t. the time step." << endl;
#ifdef P3M_TEST

    Track::block->bunch->runTests();

#else
    itsTracker = new ParallelTTracker(*Track::block->use->fetchLine(),
                                      Track::block->bunch,
                                      *ds,
                                      Track::block->reference,
                                      false,
                                      false,
                                      Track::block->localTimeSteps,
                                      Track::block->zstart,
                                      Track::block->zstop,
                                      Track::block->dT);
#endif
// #endif
    // itsTracker->setMpacflg(mpacflg); // set multipacting flag in ParallelTTracker
}

void TrackRun::setupCyclotronTracker(){
    OpalData::getInstance()->setInOPALCyclMode();
    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            BoundaryGeometry *bg = BoundaryGeometry::find(Attributes::getString(itsAttr[BOUNDARYGEOMETRY]))->
                clone(getOpalName() + std::string("_geometry"));
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }

    setupFieldsolver();

    Track::block->bunch->PType = ParticleType::REGULAR;

    std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    if (distr_str.size() == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distr_str.at(0));
    }

    // set macromass and charge for simulation particles
    double macromass = 0.0;
    double macrocharge = 0.0;

    const int specifiedNumBunch = int(std::abs(Round(Attributes::getReal(itsAttr[TURNS]))));

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (opal->hasBunchAllocated()) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    if(beam->getNumberOfParticles() < 3 || beam->getCurrent() == 0.0) {
        macrocharge = beam->getCharge() * q_e;
        macromass = beam->getMass();
        dist->createOpalCycl(Track::block->bunch,
                             beam->getNumberOfParticles(),
                             beam->getCurrent(),*Track::block->use->fetchLine());

    } else {

        /**
           getFrequency() gets RF frequency [MHz], NOT isochronous  revolution frequency of particle!
           getCurrent() gets beamcurrent [A]

        */
        macrocharge = beam->getCurrent() / (beam->getFrequency() * 1.0e6); // [MHz]-->[Hz]

        if(!opal->hasBunchAllocated()) {
            if(!opal->inRestartRun()) {
                macrocharge /= beam->getNumberOfParticles();
                macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
                dist->createOpalCycl(Track::block->bunch,
                                     beam->getNumberOfParticles(),
                                     beam->getCurrent(),
                                     *Track::block->use->fetchLine());

            } else {
                dist->doRestartOpalCycl(Track::block->bunch,
                                        beam->getNumberOfParticles(),
                                        opal->getRestartStep(),
                                        specifiedNumBunch,
                                        phaseSpaceSink_m);
                macrocharge /= beam->getNumberOfParticles();
                macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
            }
        }
    }
    Track::block->bunch->setMass(macromass); // set the Mass per macro-particle, [GeV/c^2]
    Track::block->bunch->setCharge(macrocharge);  // set the charge per macro-particle, [C]

    *gmsg << "* Mass of simulation particle= " << macromass << " GeV/c^2" << endl;
    *gmsg << "* Charge of simulation particle= " << macrocharge << " [C]" << endl;

    Track::block->bunch->setdT(1.0 / (Track::block->stepsPerTurn * beam->getFrequency() * 1.0e6));
    Track::block->bunch->setStepsPerTurn(Track::block->stepsPerTurn);

    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    if(!opal->inRestartRun())
        if(!opal->hasDataSinkAllocated()) {
            ds = new DataSink(phaseSpaceSink_m);
            opal->setDataSink(ds);
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    else {
        ds = new DataSink(phaseSpaceSink_m, -1);
        opal->setDataSink(ds);
    }

    if(!opal->hasBunchAllocated()) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == CYCLOTRON-T, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == CYCLOTRON-T, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }
    *gmsg << "* Number of neighbour bunches= " << specifiedNumBunch << endl;
    *gmsg << "* DT                         = " << Track::block->dT.front() << endl;
    *gmsg << "* MAXSTEPS                   = " << Track::block->localTimeSteps.front() << endl;
    *gmsg << "* Phase space dump frequency = " << Options::psDumpFreq << endl;
    *gmsg << "* Statistics dump frequency  = " << Options::statDumpFreq << " w.r.t. the time step." << endl;
    *gmsg << "* ********************************************************************************** " << endl;
    
    ds->setMaxNumBunches(specifiedNumBunch);

    itsTracker = new ParallelCyclotronTracker(*Track::block->use->fetchLine(),
                                              Track::block->bunch, *ds, Track::block->reference,
                                              false, false, Track::block->localTimeSteps.front(),
                                              Track::block->timeIntegrator, specifiedNumBunch);

    if(opal->inRestartRun()) {

        H5PartWrapperForPC *h5pw = static_cast<H5PartWrapperForPC*>(phaseSpaceSink_m);
        itsTracker->setBeGa(h5pw->getMeanMomentum());

        itsTracker->setPr(h5pw->getReferencePr());
        itsTracker->setPt(h5pw->getReferencePt());
        itsTracker->setPz(h5pw->getReferencePz());

        itsTracker->setR(h5pw->getReferenceR());
        itsTracker->setTheta(h5pw->getReferenceT());
        itsTracker->setZ(h5pw->getReferenceZ());

        // The following is for restarts in local frame
        itsTracker->setPhi(h5pw->getAzimuth());
        itsTracker->setPsi(h5pw->getElevation());
        itsTracker->setPreviousH5Local(h5pw->getPreviousH5Local());
    }

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    // *gmsg << *Track::block->bunch  << endl;

    if(specifiedNumBunch > 1) {

        // only for regular  run of multi bunches, instantiate the  PartBins class
        // note that for restart run of multi bunches, PartBins class is instantiated in function
        // Distribution::doRestartOpalCycl()
        if(!opal->inRestartRun()) {

            // already exist bins number initially
            const int BinCount = 1;
            //allowed maximal bin
            const int MaxBinNum = 1000;

            // initialize particles number for each bin (both existed and not yet emmitted)
            size_t partInBin[MaxBinNum];
            for(int ii = 0; ii < MaxBinNum; ii++) partInBin[ii] = 0;
            partInBin[0] =  beam->getNumberOfParticles();

            Track::block->bunch->setPBins(new PartBinsCyc(MaxBinNum, BinCount, partInBin));
            // the allowed maximal bin number is set to 100
            //Track::block->bunch->setPBins(new PartBins(100));

        }

        // mode of generating new bunches:
        // "FORCE" means generating one bunch after each revolution, until get "TURNS" bunches.
        // "AUTO" means only when the distance between two neighbor bunches is below the limitation,
        //        then starts to generate new bunches after each revolution,until get "TURNS" bunches;
        //        otherwise, run single bunch track

        *gmsg << "***---------------------------- MULTI-BUNCHES MULTI-ENERGY-BINS MODE "
              << "----------------------------*** " << endl;

        double paraMb = Attributes::getReal(itsAttr[PARAMB]);
        itsTracker->setParaAutoMode(paraMb);

        if(opal->inRestartRun()) {

            itsTracker->setLastDumpedStep(opal->getRestartStep());

            if(Track::block->bunch->pbin_m->getLastemittedBin() < 2) {
                *gmsg << "In this restart job, the multi-bunches mode is forcely set to AUTO mode." << endl;
                itsTracker->setMultiBunchMode("AUTO");
            } else {
                *gmsg << "In this restart job, the multi-bunches mode is forcely set to FORCE mode." << endl
                      << "If the existing bunch number is less than the specified number of TURN, "
                      << "readin the phase space of STEP#0 from h5 file consecutively" << endl;
                itsTracker->setMultiBunchMode("FORCE");
            }
        } else {
            std::string mbmode = Util::toUpper(Attributes::getString(itsAttr[MBMODE]));
            itsTracker->setMultiBunchMode(mbmode);
        }
        
        ParallelCyclotronTracker* tracker_p = dynamic_cast<ParallelCyclotronTracker*>(itsTracker);
        
        tracker_p->setMultiBunchEta(Attributes::getReal(itsAttr[MB_ETA]));
        
        tracker_p->setMultiBunchBinning(Attributes::getString(itsAttr[MB_BINNING]));
    }
}

void TrackRun::setupStatisticalErrors(const std::string & method) {
    std::vector<std::string> tmp;
    std::string arguments = method.substr(19,method.length() - 20);
    boost::algorithm::split(tmp, arguments, boost::algorithm::is_any_of(","));

    for (std::string &arg: tmp) {
        boost::algorithm::trim(arg);
    }

    if (Util::toUpper(tmp.at(0)) == "PARALLEL-T" ||
        Util::toUpper(tmp.at(0)) == "CYCLOTRON-T" ||
        Util::toUpper(tmp.at(0)) == "PARALLEL-SLICE" ||
        Util::toUpper(tmp.at(0)) == "THIN" ||
        Util::toUpper(tmp.at(0)) == "THICK") {

        if (tmp.size() != 3) {
            throw OpalException("TrackRun::setupStatisticalErrors()",
                                "number of arguments: " + std::to_string(tmp.size()) + " != 3");
        }

        if(!opal->hasBunchAllocated()) {
            itsTracker = new StatisticalErrors(*Track::block->use->fetchLine(),
                                               Track::block->reference,
                                               false,
                                               false,
                                               tmp.at(0),
                                               std::stoul(tmp.at(1)),
                                               std::stoul(tmp.at(2)),
                                               Attributes::getStringArray(itsAttr[OBJECTIVES]));
        } else {
            itsTracker = new NilTracker(*Track::block->use->fetchLine(),
                                        Track::block->reference,
                                        false,
                                        false);
        }

    } else {
        throw OpalException("TrackRun::setupStatisticalErrors()",
                            "unkonwn method '" + tmp.at(0) + "' provided");
    }
}

void TrackRun::setupFieldsolver() {
    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));

    if (Util::toUpper(fs->getType()) != std::string("NONE")) {
        size_t numGridPoints = fs->getMX()*fs->getMY()*fs->getMT(); // total number of gridpoints
        Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
        size_t numParticles = beam->getNumberOfParticles();

        if (!opal->inRestartRun() && numParticles < numGridPoints
#ifdef ENABLE_AMR
            && !Options::amr)
#else
            )
#endif
        {
            throw OpalException("TrackRun::setupFieldsolver()",
                                "The number of simulation particles (" + std::to_string(numParticles) + ") \n" +
                                "is smaller than the number of gridpoints (" + std::to_string(numGridPoints) + ").\n" +
                                "Please increase the number of particles or reduce the size of the mesh.\n");
        }

        OpalData::getInstance()->addProblemCharacteristicValue("MX", fs->getMX());
        OpalData::getInstance()->addProblemCharacteristicValue("MY", fs->getMY());
        OpalData::getInstance()->addProblemCharacteristicValue("MT", fs->getMT());

    }

    fs->initCartesianFields();
    Track::block->bunch->setSolver(fs);
    if (fs->hasPeriodicZ())
        Track::block->bunch->setBCForDCBeam();
    else
        Track::block->bunch->setBCAllOpen();
}

double TrackRun::setDistributionParallelT(Beam *beam) {

    // If multipacting flag is not set, get distribution(s).
    if (!Attributes::getBool(itsAttr[MULTIPACTING])) {
        /*
         * Distribution(s) can be set via a single distribution or a list
         * (array) of distributions. If an array is defined the first in the
         * list is treated as the primary distribution. All others are added to
         * it to create the full distribution.
         */
        std::vector<std::string> distributionArray
            = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
        const size_t numberOfDistributions = distributionArray.size();

        if (numberOfDistributions == 0) {
            dist = Distribution::find(defaultDistribution);
        } else {
            dist = Distribution::find(distributionArray.at(0));
            dist->setNumberOfDistributions(numberOfDistributions);

            if (numberOfDistributions > 1) {
                *gmsg << endl
                      << "---------------------------------" << endl
                      << "Found more than one distribution:" << endl << endl;
                *gmsg << "Main Distribution" << endl
                      << "---------------------------------" << endl
                      << distributionArray.at(0) << endl << endl
                      << "Secondary Distribution(s)" << endl
                      << "---------------------------------" << endl;

                for (size_t i = 1; i < numberOfDistributions; ++ i) {
                    Distribution *distribution = Distribution::find(distributionArray.at(i));
                    distribution->setNumberOfDistributions(numberOfDistributions);
                    distrs_m.push_back(distribution);

                    *gmsg << distributionArray.at(i) << endl;
                }
                *gmsg << endl
                      << "---------------------------------" << endl << endl;
            }
        }
    }

    /*
     * Initialize distributions.
     */
    size_t numberOfParticles = beam->getNumberOfParticles();
    if (!opal->hasBunchAllocated()) {
        if (!opal->inRestartRun()) {
            if (!Attributes::getBool(itsAttr[MULTIPACTING])) {
                /*
                 * Here we are not doing a restart or doing a mulitpactor run
                 * and we do not have a bunch already allocated.
                 */
                Track::block->bunch->setDistribution(dist,
                                                     distrs_m,
                                                     numberOfParticles);

                /*
                 * If this is an injected beam (rather than an emitted beam), we
                 * make sure it doesn't have any particles at z < 0.
                 */

                opal->setGlobalPhaseShift(0.5 * dist->getTEmission() + dist->getEmissionTimeShift());
            }
        } else {
            /*
             * Read in beam from restart file.
             */
            dist->doRestartOpalT(Track::block->bunch, numberOfParticles, opal->getRestartStep(), phaseSpaceSink_m);
        }
    }

    // Return charge per macroparticle.
    return beam->getCharge() * beam->getCurrent() / (beam->getFrequency()*1.0e6) / numberOfParticles;

}