// ------------------------------------------------------------------------
// $RCSfile: TrackCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackCmd
//   The class for the OPAL TRACK command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackCmd.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Track/Track.h"
#include "Track/TrackParser.h"


// Class Track
// ------------------------------------------------------------------------

namespace {

    // The attributes of class TrackRun.
    enum {
        LINE,         // The name of lattice to be tracked.
        BEAM,         // The name of beam to be used.
        DT,           // The integration timestep in second.
                      // In case of the adaptive integrator, time step guideline for
                      // external field integration.
        DTSCINIT,     // Only for adaptive integrator: Initial time step for space charge integration.
        DTAU,         // Only for adaptive integrator: Alternative way to set accuracy of space
                      // charge integration. Has no direct interpretation like DTSCINIT, but lower
                      // means smaller steps and more accurate. If given, DTSCINIT is not used. Useful
                      // for continuing with same step size in follow-up tracks.
        T0,           // The elapsed time (sec) of the bunch
        MAXSTEPS,     // The maximum timesteps we integrate
        ZSTART,       // Defines a z-location [m] where the reference particle starts
        ZSTOP,        // Defines a z-location [m], after which the simulation stops when the last particles passes
        STEPSPERTURN, // Return the timsteps per revolution period. ONLY available for OPAL-cycl.
        TIMEINTEGRATOR, // the name of time integrator
        MAP_ORDER,    // Truncation order of maps for ThickTracker (default: 1 (linear))
        SIZE
    };
}

TrackCmd::TrackCmd():
    Action(SIZE, "TRACK",
           "The \"TRACK\" command initiates tracking.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of lattice to be tracked");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "Name of beam to be used", "UNNAMED_BEAM");
    itsAttr[DT] = Attributes::makeRealArray
                  ("DT", "THE INTEGRATION TIMESTEP IN SECONDS");
    itsAttr[DTSCINIT] = Attributes::makeReal
                  ("DTSCINIT", "Only for adaptive integrator: Initial time step for space charge integration", 1e-12);
    itsAttr[DTAU] = Attributes::makeReal
                  ("DTAU", "Only for adaptive integrator: Alternative way to set accuracy of space integration.", -1.0);
    itsAttr[T0] = Attributes::makeReal
                  ("T0", "THE ELAPSED TIME OF THE BUNCH IN SECONDS", 0.0);
    itsAttr[MAXSTEPS] = Attributes::makeRealArray
                        ("MAXSTEPS", "THE MAXIMUM NUMBER OF INTEGRATION STEPS DT, should be larger ZSTOP/(beta*c average)");
    itsAttr[STEPSPERTURN] = Attributes::makeReal
                            ("STEPSPERTURN", "THE TIME STEPS PER REVOLUTION PERIOD, ONLY FOR OPAL-CYCL", 720);
    itsAttr[ZSTART] = Attributes::makeReal
                      ("ZSTART", "Defines a z-location [m] where the reference particle starts", 0.0);
    itsAttr[ZSTOP] = Attributes::makeRealArray
                     ("ZSTOP", "Defines a z-location [m], after which the simulation stops when the last particles passes");
    itsAttr[TIMEINTEGRATOR] = Attributes::makeString
                              ("TIMEINTEGRATOR", "Name of time integrator to be used", "RK-4");
    
    itsAttr[MAP_ORDER] = Attributes::makeReal
                     ("MAP_ORDER", "Truncation order of maps for ThickTracker (default: 1, i.e. linear)", 1);

    registerOwnership(AttributeHandler::COMMAND);
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "RUN");
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "START");
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "TSAVE");
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "ENDTRACK");
}

TrackCmd::TrackCmd(const std::string &name, TrackCmd *parent):
    Action(name, parent)
{}


TrackCmd::~TrackCmd()
{}


TrackCmd *TrackCmd::clone(const std::string &name) {
    return new TrackCmd(name, this);
}

std::vector<double> TrackCmd::getDT() const {
    std::vector<double> dt = Attributes::getRealArray(itsAttr[DT]);
    if (dt.size() == 0) {
        dt.push_back(1e-12);
    }
    return dt;
}

double TrackCmd::getDTSCINIT() const {
    return Attributes::getReal(itsAttr[DTSCINIT]);
}

double TrackCmd::getDTAU() const {
    return Attributes::getReal(itsAttr[DTAU]);
}

double TrackCmd::getT0() const {
    return Attributes::getReal(itsAttr[T0]);
}

double TrackCmd::getZSTART() const {
    double zstart = Attributes::getReal(itsAttr[ZSTART]);
    return zstart;
}

std::vector<double> TrackCmd::getZSTOP() const {
    std::vector<double> zstop = Attributes::getRealArray(itsAttr[ZSTOP]);
    if (zstop.size() == 0) {
        zstop.push_back(1000000.0);
    }
    return zstop;
}

std::vector<unsigned long long> TrackCmd::getMAXSTEPS() const {
    std::vector<double> maxsteps_d = Attributes::getRealArray(itsAttr[MAXSTEPS]);
    std::vector<unsigned long long> maxsteps_i;
    if (maxsteps_d.size() == 0) {
        maxsteps_i.push_back(10ul);
    }
    for (auto it = maxsteps_d.begin(); it != maxsteps_d.end(); ++ it) {
        if (*it < 0) {
            maxsteps_i.push_back(10);
        } else {
            unsigned long long value = *it;
            maxsteps_i.push_back(value);
        }
    }

    return maxsteps_i;
}

int TrackCmd::getSTEPSPERTURN() const {
    return (int) Attributes::getReal(itsAttr[STEPSPERTURN]);
}

// return int type rathor than string to improve the speed
int TrackCmd::getTIMEINTEGRATOR() const {
    std::string name = Attributes::getString(itsAttr[TIMEINTEGRATOR]);
    int  nameID;
    if(name == std::string("RK-4"))
        nameID =  0;
    else if(name == std::string("LF-2"))
        nameID =  1;
    else if(name == std::string("MTS"))
        nameID = 2;
    else if(name == std::string("AMTS"))
        nameID = 3;
    else
        nameID = -1;

    return nameID;
}

void TrackCmd::execute() {
    // Find BeamSequence and Beam definitions.
    BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    std::vector<double> dt = getDT();
    double t0 = getT0();
    std::vector<unsigned long long> maxsteps = getMAXSTEPS();
    int    stepsperturn = getSTEPSPERTURN();
    double zstart = getZSTART();
    std::vector<double> zstop = getZSTOP();
    int timeintegrator = getTIMEINTEGRATOR();

    size_t numTracks = dt.size();
    numTracks = std::max(numTracks, maxsteps.size());
    numTracks = std::max(numTracks, zstop.size());
    for (size_t i = dt.size(); i < numTracks; ++ i) {
        dt.push_back(dt.back());
    }
    for (size_t i = maxsteps.size(); i < numTracks; ++ i) {
        maxsteps.push_back(maxsteps.back());
    }
    for (size_t i = zstop.size(); i < numTracks; ++ i) {
        zstop.push_back(zstop.back());
    }

   // Execute track block.
    Track::block = new Track(use, beam->getReference(), dt, maxsteps,
                             stepsperturn, zstart, zstop,
                             timeintegrator, t0, getDTSCINIT(), getDTAU());
    
    Track::block->truncOrder = (int)Attributes::getReal(itsAttr[MAP_ORDER]);
    
    Track::block->parser.run();

    // Clean up.
    delete Track::block;
    Track::block = 0;
}