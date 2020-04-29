//
// Class Track
//   Hold data for tracking.
//   Acts as a communication area between the various tracking commands.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#ifndef OPAL_Track_HH
#define OPAL_Track_HH

#include "Algorithms/PartData.h"
#include "Track/TrackParser.h"
#include <stack>
#include <vector>

class BeamSequence;

template <class T, unsigned Dim>
class PartBunchBase;

class Track {

public:

    Track(BeamSequence *, const PartData &, const std::vector<double> & dt,
          const std::vector<unsigned long long> & maxtsteps, int stepsperturn,
          double zStart, const std::vector<double> & zStop, int timeintegrator,
          double t0, double dtScInit, double deltaTau);
    ~Track();

    /// The particle bunch to be tracked.
    PartBunchBase<double, 3> *bunch;

    /// The reference data.
    PartData reference;

    /// The lattice to be tracked through.
    BeamSequence *use;

    /// The parser used during tracking.
    TrackParser parser;

    /// The block of track data.
    static Track *block;

    static void stash();
    static Track* pop();

    /// The initial timestep
    std::vector<double> dT;

    // For AMTS integrator in OPAL-T
    double dtScInit, deltaTau;

    /// The ellapsed time of the beam can be used to propper
    /// start the beam when created in a cavity i.e. without emission
    double t0_m;

    /// Maximal number of timesteps
    std::vector<unsigned long long> localTimeSteps;

    /// The timsteps per revolution period. ONLY available for OPAL-cycl.
    int stepsPerTurn;

    /// The location at which the simulation starts
    double zstart;

    /// The location at which the simulation stops
    std::vector<double> zstop;

    /// The ID of time integrator
    // 0 --- RK-4(default)
    // 1 --- LF-2
    // 2 --- MTS
    // 3 --- AMTS
    int timeIntegrator;
    
    /// Trunction order for map tracking
    int truncOrder;
    
private:

    // Not implemented.
    Track();
    Track(const Track &);
    void operator=(const Track &);

    static std::stack<Track*> stashedTrack;
};

#endif // OPAL_Track_HH