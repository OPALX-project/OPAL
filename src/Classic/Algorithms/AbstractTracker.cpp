// ------------------------------------------------------------------------
// $RCSfile: AbstractTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractTracker
//   An abstract visitor class allowing to tracking a bunch of particles
//   through a beamline.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractTracker.h"
#include "Algorithms/PartData.h"


//: Abstract tracker class.
//  An abstract visitor class implementing the default behaviour for all
//  visitors capable of tracking a particle bunch through a beam line.
//  It implements access to the bunch, and keeps track of the beam
//  reference data.

// Class AbstractTracker
// ------------------------------------------------------------------------


AbstractTracker::AbstractTracker(const Beamline &beamline,
                                 const PartData &reference,
                                 bool backBeam, bool backTrack):
    DefaultVisitor(beamline, backBeam, backTrack),
    itsReference(reference)
{}


AbstractTracker::~AbstractTracker()
{}