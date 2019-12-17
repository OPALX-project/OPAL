#ifndef CLASSIC_BeamlineVisitor_HH
#define CLASSIC_BeamlineVisitor_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamlineVisitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamlineVisitor
//   Defines the abstract interface for a BeamlineVisitor.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

// Generic element classes interacting with a BeamlineVisitor.
class Component;

// Beam line structure classes.
class Beamline;
class AlignWrapper;
class CorrectorWrapper;
class FlaggedElmPtr;
class MultipoleWrapper;
class RBendWrapper;
class SBendWrapper;
class CyclotronWrapper;

// Specific element classes interacting with a BeamlineVisitor
class BeamBeam;
class BeamStripping;
class CCollimator;
class Corrector;
class Cyclotron;
class CyclotronValley;
class Degrader;
class Diagnostic;
class Drift;
class Undulator;
class FlexibleCollimator;
class Lambertson;
class Marker;
class Monitor;
class Multipole;
class MultipoleT;
class MultipoleTStraight;
class MultipoleTCurvedConstRadius;
class MultipoleTCurvedVarRadius;
class Offset;
class ParallelPlate;
class Patch;
class Probe;
class RBend;
class RBend3D;
class RFCavity;
class RFQuadrupole;
class Ring;
class SBend;
class SBend3D;
class ScalingFFAMagnet;
class Separator;
class Septum;
class Solenoid;
class Source;
class Stripper;
class TravelingWave;
class VariableRFCavity;
class VariableRFCavityFringeField;
class VerticalFFAMagnet;

// Integrators.
class Integrator;
class MapIntegrator;
class TrackIntegrator;


// Class BeamlineVisitor
// ------------------------------------------------------------------------
/// Abstract algorithm.
//  The abstract class BeamlineVisitor is the base class for all visitors
//  (algorithms) that can iterator over a beam line representation.
//  A BeamlineVisitor applies itself to the representation via the
//  ``Visitor'' pattern, see
//  [p]
//  E. Gamma, R. Helm, R. Johnson, and J. Vlissides,
//  [BR]
//  Design Patterns, Elements of Reusable Object-Oriented Software.
//  [p]
//  By using only pure abstract classes as an interface between the
//  BeamlineVisitor and the beam line representation,
//  we decouple the former from the implementation details of the latter.
//  [p]
//  The interface is defined in such a way that a visitor cannot modify the
//  structure of a beam line, but it can assign special data like misalignments
//  or integrators without problems.

class BeamlineVisitor {

public:

    BeamlineVisitor();
    virtual ~BeamlineVisitor();

    /// Execute the algorithm on the attached beam line.
    virtual void execute() = 0;

    /// Apply the algorithm to a beam-beam interaction.
    virtual void visitBeamBeam(const BeamBeam &) = 0;

    /// Apply the algorithm to a beam stripping interaction.
    virtual void visitBeamStripping(const BeamStripping &) = 0;

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &) = 0;

    /// Apply the algorithm to an arbitrary component (catch all).
    virtual void visitComponent(const Component &) = 0;

    /// Apply the algorithm to a closed orbit corrector.
    virtual void visitCorrector(const Corrector &) = 0;

    /// Apply the algorithm to a diagnostic.
    virtual void visitDegrader(const Degrader &) = 0;

    /// Apply the algorithm to a diagnostic.
    virtual void visitDiagnostic(const Diagnostic &) = 0;

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift &) = 0;

    /// Apply the algorithm to a undulator space.
    virtual void visitUndulator(const Undulator &) = 0;

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &) = 0;

    /// Apply the algorithm to an Ring
    virtual void visitRing(const Ring &) = 0;

    /// Apply the algorithm to a cyclotron.
    virtual void visitCyclotron(const Cyclotron &) = 0;

    /// Apply the algorithm to a Lambertson septum magnet.
    virtual void visitLambertson(const Lambertson &) = 0;

    /// Apply the algorithm to an Offset (placement).
    virtual void visitOffset(const Offset &) = 0;

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &) = 0;

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor &) = 0;

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &) = 0;

    /// Apply the algorithm to an arbitrary Multipole.
    virtual void visitMultipoleT(const MultipoleT &) = 0;

    /// Apply the algorithm to an arbitrary straight Multipole.
    virtual void visitMultipoleTStraight(const MultipoleTStraight &) = 0;

    /// Apply the algorithm to an arbitrary curved Multipole of constant radius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &) = 0;

    /// Apply the algorithm to an arbitrary curved Multipole of variable radius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &) = 0;

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &) = 0;

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &) = 0;

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &) = 0;

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &) = 0;

    /// Apply the algorithm to a variable RF cavity.
    virtual void visitVariableRFCavity(const VariableRFCavity &) = 0;

    /// Apply the algorithm to a variable RF cavity with Fringe Field.
    virtual void visitVariableRFCavityFringeField
                                      (const VariableRFCavityFringeField &) = 0;

    /// Apply the algorithm to a RF cavity.
    virtual void visitTravelingWave(const TravelingWave &) = 0;

    /// Apply the algorithm to a RF quadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &) = 0;

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &) = 0;

    /// Apply the algorithm to a Sector Bend with 3D field map.
    virtual void visitSBend3D(const SBend3D &) = 0;

    /// Apply the algorithm to an electrostatic separator.
    virtual void visitSeparator(const Separator &) = 0;

    /// Apply the algorithm to a septum magnet.
    virtual void visitSeptum(const Septum &) = 0;

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &) = 0;

    /// Apply the algorithm to a solenoid.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet &) = 0;

    /// Apply the algorithm to a source.
    virtual void visitSource(const Source &) = 0;

    /// Apply the algorithm to a Beamline.
    virtual void visitBeamline(const Beamline &) = 0;

    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &) = 0;

    /// Apply the algorithm to an align wrapper.
    virtual void visitAlignWrapper(const AlignWrapper &) = 0;

    /// Apply the algorithm to an corrector wrapper.
    virtual void visitCorrectorWrapper(const CorrectorWrapper &) = 0;

    /// Apply the algorithm to an corrector wrapper.
    virtual void visitCyclotronWrapper(const CyclotronWrapper &) = 0;

    /// Apply the algorithm to an multipole wrapper.
    virtual void visitMultipoleWrapper(const MultipoleWrapper &) = 0;

    /// Apply the algorithm to an RBend wrapper.
    virtual void visitRBendWrapper(const RBendWrapper &) = 0;

    /// Apply the algorithm to an SBend wrapper.
    virtual void visitSBendWrapper(const SBendWrapper &) = 0;

    /// Apply the algorithm to an ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &) = 0;

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &) = 0;

    /// Apply the algorithm to a particle stripper.
    virtual void visitStripper(const Stripper &) = 0;

    /// Apply the algorithm to a generic integrator.
    virtual void visitIntegrator(const Integrator &) = 0;

    /// Apply the algorithm to an integrator capable of tracking.
    virtual void visitTrackIntegrator(const TrackIntegrator &) = 0;

    /// Apply the algorithm to an integrator capable of mapping.
    virtual void visitMapIntegrator(const MapIntegrator &) = 0;

    /// Apply the algorithm to a vertical FFA magnet
    virtual void visitVerticalFFAMagnet(const VerticalFFAMagnet &) = 0;

private:

    // Not implemented.
    BeamlineVisitor(const BeamlineVisitor &);
    void operator=(const BeamlineVisitor &);
};

inline
void BeamlineVisitor::visitRBend3D(const RBend3D &) {

}

#endif // CLASSIC_BeamlineVisitor_HH
