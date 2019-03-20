#ifndef CLASSIC_DefaultVisitor_HH
#define CLASSIC_DefaultVisitor_HH

// ------------------------------------------------------------------------
// $RCSfile: DefaultVisitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DefaultVisitor
//   The default interface for a BeamlineVisitor.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:04 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamlineVisitor.h"

class ElementBase;

// Class DefaultVisitor
// ------------------------------------------------------------------------
/// Default algorithms.
//  A default implementation for all visitors that can iterate over a
//  beam line representation.
//  This abstract base class implements the default behaviour for the
//  structural classes Beamline and FlaggedElmPtr, and for all wrappers.
//  It also holds the data required for all visitors in a protected area.

class DefaultVisitor: public BeamlineVisitor {

public:

    /// Constructor.
    //  Arguments:
    //  [ol]
    //  [li]The beamline to be used.
    //  [li]If true, the beam runs backwards through the line.
    //  [li]If true, we track against the beam.
    //  [/ol]
    DefaultVisitor(const Beamline &beamline, bool backBeam, bool backTrack);

    virtual ~DefaultVisitor() = 0;

    /// Apply the algorithm to the top-level beamline.
    virtual void execute();

    /// Apply the algorithm to a beam-beam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to an cyclotron
    virtual void visitCyclotron(const Cyclotron &);

    /// Apply the algorithm to an opal ring..
    virtual void visitRing(const Ring &);

    /// Apply the algorithm to a corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a drift.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a drift.
    virtual void visitUndulator(const Undulator &);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a Lambertson.
    virtual void visitLambertson(const Lambertson &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to a multipoleT.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to a multipoleTStraight.
    virtual void visitMultipoleTStraight(const MultipoleTStraight &);

    /// Apply the algorithm to a multipoleTCurvedConstRadius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &);

    /// Apply the algorithm to a multipoleTCurvedVarRadius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &);

    /// Apply the algorithm to an Offset.
    virtual void visitOffset(const Offset &);

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitVariableRFCavity(const VariableRFCavity &vcav);

    /// Apply the algorithm to a RF cavity.
    virtual void visitVariableRFCavityFringeField(const VariableRFCavityFringeField &vcav);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitTravelingWave(const TravelingWave &);

    /// Apply the algorithm to a RF quadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend3D(const SBend3D &);

    /// Apply the algorithm to a scaling FFA magnet.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet &);

    /// Apply the algorithm to a separator.
    virtual void visitSeparator(const Separator &);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a source.
    virtual void visitSource(const Source &);

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);

    /// Apply the algorithm to a charge stripper.
    virtual void visitStripper(const Stripper &);

    /// Apply the algorithm to a beam line.
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);


    /// Apply the algorithm to an align wrapper..
    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to an corrector wrapper..
    virtual void visitCorrectorWrapper(const CorrectorWrapper &);

    /// Apply the algorithm to an cyclotron wrapper..
    virtual void visitCyclotronWrapper(const CyclotronWrapper &);

    /// Apply the algorithm to an multipole wrapper..
    virtual void visitMultipoleWrapper(const MultipoleWrapper &);

    /// Apply the algorithm to an RBend wrapper..
    virtual void visitRBendWrapper(const RBendWrapper &);

    /// Apply the algorithm to an SBend wrapper..
    virtual void visitSBendWrapper(const SBendWrapper &);


    /// Apply the algorithm to a generic integrator.
    virtual void visitIntegrator(const Integrator &);

    /// Apply the algorithm to an integrator capable of mapping.
    virtual void visitMapIntegrator(const MapIntegrator &);

    /// Apply the algorithm to an integrator capable of tracking.
    virtual void visitTrackIntegrator(const TrackIntegrator &);

protected:

    // The top level beamline.
    const Beamline &itsLine;

    // The direction flags and corresponding factors.
    bool back_beam;   // true, if beam runs from right (s=C) to left (s=0).
    bool back_track;  // true, if tracking opposite to the beam direction.
    bool back_path;   // true, if tracking from right (s=C) to left (s=0).
    // back_path = back_beam && ! back_track || back_track && ! back_beam.

    double flip_B;    // set to -1.0 to flip B fields, when back_beam is true.
    double flip_s;    // set to -1.0 to flip direction of s,
    // when back_path is true.

private:

    // Not implemented.
    DefaultVisitor();
    DefaultVisitor(const DefaultVisitor &);
    void operator=(const DefaultVisitor &);

    // Default do-nothing routine.
    virtual void applyDefault(const ElementBase &);

    // The element order flag. Initially set to back_path.
    // This flag is reversed locally for reflected beam lines.
    bool local_flip;
};

#endif // CLASSIC_DefaultVisitor_HH
