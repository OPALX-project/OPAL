//
// Class NilTracker
//   :FIXME: Add class description
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Algorithms/NilTracker.h"

NilTracker::NilTracker(const Beamline &beamline,
                       const PartData &reference,
                       bool revBeam,
                       bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack)
{ }


NilTracker::~NilTracker() {

}

void NilTracker::execute() {

}