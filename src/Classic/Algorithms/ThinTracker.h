#ifndef CLASSIC_ThinTracker_HH
#define CLASSIC_ThinTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ThinTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThinTracker
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"

class BMultipoleField;

template <class T, unsigned Dim>
class PartBunchBase;
class PlanarArcGeometry;


// Class ThinTracker
// ------------------------------------------------------------------------
/// Track with thin lens algorithm.
//  The visitor class for tracking a bunch of particles through a beamline
//  using a thin-lens approximation for all elements.
//  [P]
//  Approximations used:
//  [UL]
//  [LI]All active elements are represented as thin lenses, sandwiched
//    between two drifts, each half of the element length.
//  [LI]Drifts are handled with a second-order approximation.
//  [LI]Geometric transformations ignore rotations about transverse axes
//    and translations along the design orbit and truncate after second order.
//  [/UL]

class ThinTracker: public Tracker {

public:

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    ThinTracker(const Beamline &bl, const PartData &data,
                bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    ThinTracker(const Beamline &bl, PartBunchBase<double, 3> *bunch,
                const PartData &data, bool revBeam, bool revTrack);

    virtual ~ThinTracker();

    /// Apply algorithm to BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply algorithm to Collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply algorithm to Corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a drift.
    virtual void visitDegrader(const Degrader &);

    /// Apply algorithm to Diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply algorithm to Drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a undulator.
    virtual void visitUndulator(const Undulator &);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply algorithm to Lambertson.
    virtual void visitLambertson(const Lambertson &);

    /// Apply algorithm to Marker.
    virtual void visitMarker(const Marker &);

    /// Apply algorithm to Monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply algorithm to Multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply algorithm to Probe.
    virtual void visitProbe(const Probe &);

    /// Apply algorithm to RBend.
    virtual void visitRBend(const RBend &);

    /// Apply algorithm to RFCavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply algorithm to RFQuadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply algorithm to SBend.
    virtual void visitSBend(const SBend &);

    /// Apply algorithm to Separator.
    virtual void visitSeparator(const Separator &);

    /// Apply algorithm to Septum.
    virtual void visitSeptum(const Septum &);

    /// Apply algorithm to Solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply algorithm to Solenoid.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);

private:

    // Not implemented.
    ThinTracker();
    ThinTracker(const ThinTracker &);
    void operator=(const ThinTracker &);

    // Apply a drift length.
    // Approximate method to gain speed.
    void applyDrift(double length);
};

#endif // CLASSIC_ThinTracker_HH
