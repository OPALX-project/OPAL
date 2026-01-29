#ifndef CLASSIC_AbstractTracker_HH
#define CLASSIC_AbstractTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: AbstractTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractTracker
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/PartData.h"

// Class AbstractTracker
// ------------------------------------------------------------------------
/// Track particles or bunches.
//  An abstract base class for all visitors capable of tracking particles
//  through a beam element.
//  This class redefines all visitXXX() methods for elements as pure
//  to force their implementation in derived classes.

class AbstractTracker: public DefaultVisitor {

public:

    // Particle coordinate numbers.
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    AbstractTracker(const Beamline &, const PartData &,
                    bool backBeam, bool backTrack);

    virtual ~AbstractTracker();

protected:

    /// The reference information.
    const PartData itsReference;

private:

    // Not implemented.
    AbstractTracker();
    AbstractTracker(const AbstractTracker &);
    void operator=(const AbstractTracker &);
};

#endif // CLASSIC_AbstractTracker_HH