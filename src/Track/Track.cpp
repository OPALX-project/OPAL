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

#include "Track/Track.h"
// #include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartBunch.h" //FIXME
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"

Track *Track::block = 0;
std::stack<Track*> Track::stashedTrack;

/**
Track is asking the dictionary if already a
particle bunch was allocated. If that is the
case Track is using the already allocated bunch,
otherwise a new bunch is allocated in the dictionary.
*/


Track::Track(BeamSequence *u, const PartData &ref, const std::vector<double> & dt,
             const std::vector<unsigned long long> & maxtsteps, int stepsperturn,
             double zStart, const std::vector<double> & zStop, int timeintegrator,
             double t0, double dtScInit, double deltaTau):
    bunch(nullptr),
    reference(ref),
    use(u),
    parser(),
    dT(dt),
    dtScInit(dtScInit),
    deltaTau(deltaTau),
    t0_m(t0),
    localTimeSteps(maxtsteps),
    stepsPerTurn(stepsperturn),
    zstart(zStart),
    zstop(zStop),
    timeIntegrator(timeintegrator),
    truncOrder(1)
{
    if(!OpalData::getInstance()->hasBunchAllocated()) {
#ifdef ENABLE_AMR
        if (Options::amr)
            OpalData::getInstance()->setPartBunch(new AmrPartBunch(&ref));
        else
#endif
            OpalData::getInstance()->setPartBunch(new PartBunch(&ref));
    }

    bunch = OpalData::getInstance()->getPartBunch();
}


Track::~Track()
{}

void Track::stash() {
    PAssert_EQ(stashedTrack.size(), 0);

    stashedTrack.push(block);
    block = 0;
}

Track* Track::pop() {
    delete block;
    block = stashedTrack.top();
    stashedTrack.pop();

    return block;
}
