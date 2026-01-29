//
// Class DefaultVisitor
//   The default interface for a BeamlineVisitor.
//   A default implementation for all visitors that can iterate over a
//   beam line representation.
//   This abstract base class implements the default behaviour for the
//   structural classes Beamline and FlaggedElmPtr.
//   It also holds the data required for all visitors in a protected area.
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
#ifndef CLASSIC_DefaultVisitor_HH
#define CLASSIC_DefaultVisitor_HH

#include "AbsBeamline/BeamlineVisitor.h"

class ElementBase;

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

    /// Apply the algorithm to a beam line.
    virtual void visitBeamline(const Beamline &);
    
    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a closed orbit corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to an cyclotron
    virtual void visitCyclotron(const Cyclotron &);

    /// Apply the algorithm to a degrader.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);
    
    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to an arbitrary straight multipole.
    virtual void visitMultipoleTStraight(const MultipoleTStraight &);

    /// Apply the algorithm to an arbitrary curved multipole of constant radius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &);

    /// Apply the algorithm to an arbitrary curved multipole of variable radius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &);

    /// Apply the algorithm to an offset (placement).
    virtual void visitOffset(const Offset &);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a ring.
    virtual void visitRing(const Ring &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a sector bend with 3D field map.
    virtual void visitSBend3D(const SBend3D &);

    /// Apply the algorithm to a scaling FFA magnet.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet &);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a source.
    virtual void visitSource(const Source &);

    /// Apply the algorithm to a particle stripper.
    virtual void visitStripper(const Stripper &);

    /// Apply the algorithm to a traveling wave.
    virtual void visitTravelingWave(const TravelingWave &);

#ifdef ENABLE_OPAL_FEL
    /// Apply the algorithm to an undulator.
    virtual void visitUndulator(const Undulator &);
#endif

    /// Apply the algorithm to a vacuum space.
    virtual void visitVacuum(const Vacuum &);

    /// Apply the algorithm to a a variable RF cavity.
    virtual void visitVariableRFCavity(const VariableRFCavity &vcav);

    /// Apply the algorithm to a a variable RF cavity with Fringe Field.
    virtual void visitVariableRFCavityFringeField(const VariableRFCavityFringeField &vcav);

    /// Apply the algorithm to a vertical FFA magnet.
    virtual void visitVerticalFFAMagnet(const VerticalFFAMagnet &);

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
