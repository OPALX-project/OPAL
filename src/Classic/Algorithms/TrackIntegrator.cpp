// ------------------------------------------------------------------------
// $RCSfile: TrackIntegrator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackIntegrator
//   An abstract base class for integrators.
//   A TrackIntegrator propagates a single particle or particle bunch.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/TrackIntegrator.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/LogicalError.h"


// Class TrackIntegrator
// ------------------------------------------------------------------------

TrackIntegrator::TrackIntegrator(ElementBase *elem):
    Integrator(elem)
{}


TrackIntegrator::TrackIntegrator(const TrackIntegrator &rhs):
    Integrator(rhs)
{}


TrackIntegrator::~TrackIntegrator()
{}


void TrackIntegrator::accept(BeamlineVisitor &visitor) const {
    visitor.visitTrackIntegrator(*this);
}


void TrackIntegrator::trackMap(FVps<double, 6> &, const PartData &,
                               bool, bool) const {
    throw LogicalError("TrackIntegrator::trackMap()",
                       "You cannot track a map using a track integrator.");
}