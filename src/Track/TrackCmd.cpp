//
// Class TrackCmd
//   The class for the OPAL TRACK command.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Track/TrackCmd.h"

#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Track/Track.h"
#include "Track/TrackParser.h"
#include "Utilities/OpalException.h"


namespace {
    // The attributes of class TrackCmd
    enum {
        LINE,           // The name of lattice to be tracked.
        BEAM,           // The name of beam to be used.
        DT,             // The integration timestep in second.
                        // In case of the adaptive integrator, time step guideline for
                        // external field integration.
        DTSCINIT,       // Only for adaptive integrator: Initial time step for space charge integration.
        DTAU,           // Only for adaptive integrator: Alternative way to set accuracy of space
                        // charge integration. Has no direct interpretation like DTSCINIT, but lower
                        // means smaller steps and more accurate. If given, DTSCINIT is not used. Useful
                        // for continuing with same step size in follow-up tracks.
        T0,             // The elapsed time (sec) of the bunch
        MAXSTEPS,       // The maximum timesteps we integrate
        ZSTART,         // Defines a z-location [m] where the reference particle starts
        ZSTOP,          // Defines a z-location [m], after which the simulation stops when the last particles passes
        STEPSPERTURN,   // Return the timsteps per revolution period. ONLY available for OPAL-cycl.
        TIMEINTEGRATOR, // the name of time integrator
        MAP_ORDER,      // Truncation order of maps for ThickTracker (default: 1 (linear))
        SIZE
    };
}

const std::map<std::string, Steppers::TimeIntegrator> TrackCmd::stringTimeIntegrator_s = {
    {"RK-4", Steppers::TimeIntegrator::RK4},
    {"RK4",  Steppers::TimeIntegrator::RK4},
    {"LF-2", Steppers::TimeIntegrator::LF2},
    {"LF2",  Steppers::TimeIntegrator::LF2},
    {"MTS",  Steppers::TimeIntegrator::MTS}
};


TrackCmd::TrackCmd():
    Action(SIZE, "TRACK",
           "The \"TRACK\" command initiates tracking.") {
    itsAttr[LINE] = Attributes::makeString
        ("LINE", "Name of lattice to be tracked.");

    itsAttr[BEAM] = Attributes::makeString
        ("BEAM", "Name of beam to be used.", "UNNAMED_BEAM");

    itsAttr[DT] = Attributes::makeRealArray
        ("DT", "The integration timestep in [s].");

    itsAttr[DTSCINIT] = Attributes::makeReal
        ("DTSCINIT", "Only for adaptive integrator: Initial time step for space charge integration.", 1e-12);

    itsAttr[DTAU] = Attributes::makeReal
        ("DTAU", "Only for adaptive integrator: Alternative way to set accuracy of space integration.", -1.0);

    itsAttr[T0] = Attributes::makeReal
        ("T0", "The elapsed time of the bunch in seconds", 0.0);

    itsAttr[MAXSTEPS] = Attributes::makeRealArray
        ("MAXSTEPS", "The maximum number of integration steps dt, should be larger ZSTOP/(beta*c average).");

    itsAttr[ZSTART] = Attributes::makeReal
        ("ZSTART", "Defines a z-location [m] where the reference particle starts.", 0.0);

    itsAttr[ZSTOP] = Attributes::makeRealArray
        ("ZSTOP", "Defines a z-location [m], after which the simulation stops when the last particles passes.");

    itsAttr[STEPSPERTURN] = Attributes::makeReal
        ("STEPSPERTURN", "The time steps per revolution period, only for opal-cycl.", 720);

    itsAttr[TIMEINTEGRATOR] = Attributes::makePredefinedString
        ("TIMEINTEGRATOR", "Name of time integrator to be used.",
         {"RK-4", "RK4", "LF-2", "LF2", "MTS"}, "RK4");

    itsAttr[MAP_ORDER] = Attributes::makeReal
        ("MAP_ORDER", "Truncation order of maps for ThickTracker (default: 1, i.e. linear).", 1);

    registerOwnership(AttributeHandler::COMMAND);
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "RUN");
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "ENDTRACK");
}

TrackCmd::TrackCmd(const std::string& name, TrackCmd* parent):
    Action(name, parent)
{}


TrackCmd::~TrackCmd()
{}


TrackCmd* TrackCmd::clone(const std::string& name) {
    return new TrackCmd(name, this);
}

std::vector<double> TrackCmd::getDT() const {
    std::vector<double> dTs = Attributes::getRealArray(itsAttr[DT]);
    if (dTs.size() == 0) {
        dTs.push_back(1e-12);
    }
    for (double dt : dTs) {
        if (dt < 0.0) {
            throw OpalException("TrackCmd::getDT",
                                "The time steps provided with DT have to be positive");
        }
    }
    return dTs;
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

double TrackCmd::getZStart() const {
    return Attributes::getReal(itsAttr[ZSTART]);
}

std::vector<double> TrackCmd::getZStop() const {
    std::vector<double> zstop = Attributes::getRealArray(itsAttr[ZSTOP]);
    if (zstop.size() == 0) {
        zstop.push_back(1000000.0);
    }
    return zstop;
}

std::vector<unsigned long long> TrackCmd::getMaxSteps() const {
    std::vector<double> maxsteps_d = Attributes::getRealArray(itsAttr[MAXSTEPS]);
    std::vector<unsigned long long> maxsteps_i;
    if (maxsteps_d.size() == 0) {
        maxsteps_i.push_back(10ul);
    }
    for (double numSteps : maxsteps_d) {
        if (numSteps < 0) {
            throw OpalException("TrackCmd::getMAXSTEPS",
                                "The number of steps provided with MAXSTEPS has to be positive");
        } else {
            unsigned long long value = numSteps;
            maxsteps_i.push_back(value);
        }
    }

    return maxsteps_i;
}

int TrackCmd::getStepsPerTurn() const {
    return (int) Attributes::getReal(itsAttr[STEPSPERTURN]);
}

Steppers::TimeIntegrator TrackCmd::getTimeIntegrator() {
    std::string name = Attributes::getString(itsAttr[TIMEINTEGRATOR]);
    return stringTimeIntegrator_s.at(name);
}

void TrackCmd::execute() {
    // Find BeamSequence and Beam definitions.
    BeamSequence* use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
    Beam* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    std::vector<double> dt = getDT();
    double t0 = getT0();
    double dtScInit = getDTSCINIT();
    double deltaTau = getDTAU();
    std::vector<unsigned long long> maxsteps = getMaxSteps();
    int stepsperturn = getStepsPerTurn();
    double zstart = getZStart();
    std::vector<double> zstop = getZStop();
    Steppers::TimeIntegrator timeintegrator = getTimeIntegrator();

    size_t numTracks = dt.size();
    numTracks = std::max(numTracks, maxsteps.size());
    numTracks = std::max(numTracks, zstop.size());
    for (size_t i = dt.size(); i < numTracks; ++ i) {
        dt.push_back(dt.back());
    }
    for (size_t i = maxsteps.size(); i < numTracks; ++i) {
        maxsteps.push_back(maxsteps.back());
    }
    for (size_t i = zstop.size(); i < numTracks; ++i) {
        zstop.push_back(zstop.back());
    }

   // Execute track block.
    Track::block = new Track(use, beam->getReference(), dt, maxsteps,
                             stepsperturn, zstart, zstop,
                             timeintegrator, t0, dtScInit, deltaTau);

    Track::block->truncOrder = (int)Attributes::getReal(itsAttr[MAP_ORDER]);

    Track::block->parser.run();

    // Clean up.
    delete Track::block;
    Track::block = 0;
}
