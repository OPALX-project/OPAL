#ifndef CLASSIC_TrackIntegrator_HH
#define CLASSIC_TrackIntegrator_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackIntegrator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackIntegrator
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Integrator.h"

class PartData;

namespace classic {
class Particle;
}

template <class T, int N> class FVps;


// Class TrackIntegrator
// ------------------------------------------------------------------------
/// Integrate particle.
//  An abstract base class for all integrators capable of tracking particles
//  or particle bunches through a beam element.  It is assumed that this
//  integrator has no maps available.  It can therefore not track a map.

class TrackIntegrator: public Integrator {

public:

    explicit TrackIntegrator(ElementBase *);
    TrackIntegrator(const TrackIntegrator &);
    virtual ~TrackIntegrator();


    /// Apply visitor.
    virtual void accept(BeamlineVisitor &visitor) const;

    /// Make a clone.
    virtual TrackIntegrator *clone() const = 0;

    /// Track a map.
    //  The map is stored in [b]map[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    //  This version throws LogicalError, since the integrator has no map.
    virtual void trackMap(FVps<double, 6> &map, const PartData &data,
                          bool backBeam, bool backTrack) const;

private:

    // Not Implemented.
    TrackIntegrator();
    void operator=(const TrackIntegrator &);
};

#endif // CLASSIC_TrackIntegrator_HH