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
#include "Algorithms/DefaultVisitor.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/Ring.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/Source.h"
#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/Vacuum.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "AbsBeamline/VerticalFFAMagnet.h"

#ifdef ENABLE_OPAL_FEL
#include "AbsBeamline/Undulator.h"
#endif

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"



DefaultVisitor::DefaultVisitor(const Beamline &beamline,
                               bool backBeam, bool backTrack):
    itsLine(beamline), back_beam(backBeam), back_track(backTrack) {
    local_flip = back_path =
        ( back_beam && ! back_track )  ||  ( back_track && ! back_beam );
    flip_B     = back_beam ? -1.0 : 1.0;
    flip_s     = back_path ? -1.0 : 1.0;
}


DefaultVisitor::~DefaultVisitor()
{}


void DefaultVisitor::execute() {
    local_flip = ( back_beam && ! back_track )  ||  ( back_track && ! back_beam );
    itsLine.accept(*this);
}


void DefaultVisitor::visitCCollimator(const CCollimator &coll) {
    applyDefault(coll);
}

void DefaultVisitor::visitComponent(const Component &comp) {
    applyDefault(comp);
}

void DefaultVisitor::visitCorrector(const Corrector &corr) {
    applyDefault(corr);
}

void DefaultVisitor::visitCyclotron(const Cyclotron &cyc) {
    applyDefault(cyc);
}

void DefaultVisitor::visitDegrader(const Degrader &deg) {
    applyDefault(deg);
}

void DefaultVisitor::visitDrift(const Drift &drf) {
    applyDefault(drf);
}

void DefaultVisitor::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDefault(coll);
}

void DefaultVisitor::visitMarker(const Marker &mark) {
    applyDefault(mark);
}

void DefaultVisitor::visitMonitor(const Monitor &mon) {
    applyDefault(mon);
}

void DefaultVisitor::visitMultipole(const Multipole &mult) {
    applyDefault(mult);
}

void DefaultVisitor::visitMultipoleT(const MultipoleT &multT) {
    applyDefault(multT);
}

void DefaultVisitor::visitMultipoleTStraight(const MultipoleTStraight &multTstraight) {
    applyDefault(multTstraight);
}

void DefaultVisitor::visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &multTccurv) {
    applyDefault(multTccurv);
}

void DefaultVisitor::visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &multTvcurv) {
    applyDefault(multTvcurv);
}

void DefaultVisitor::visitOffset(const Offset& off) {
    applyDefault(off);
}

void DefaultVisitor::visitProbe(const Probe &probe) {
    applyDefault(probe);
}

void DefaultVisitor::visitRBend(const RBend &bend) {
    applyDefault(bend);
}

void DefaultVisitor::visitRBend3D(const RBend3D &bend) {
    applyDefault(bend);
}

void DefaultVisitor::visitRFCavity(const RFCavity &cav) {
    applyDefault(cav);
}

void DefaultVisitor::visitRing(const Ring &ring) {
   applyDefault(ring);
}

void DefaultVisitor::visitSBend(const SBend &bend) {
    applyDefault(bend);
}

void DefaultVisitor::visitSBend3D(const SBend3D &bend) {
    applyDefault(bend);
}

void DefaultVisitor::visitScalingFFAMagnet(const ScalingFFAMagnet &spiral) {
    applyDefault(spiral);
}

void DefaultVisitor::visitSeptum(const Septum &sept) {
    applyDefault(sept);
}

void DefaultVisitor::visitSolenoid(const Solenoid &sol) {
    applyDefault(sol);
}

void DefaultVisitor::visitSource(const Source &sou) {
    applyDefault(sou);
}

void DefaultVisitor::visitStripper(const Stripper &stripper) {
    applyDefault(stripper);
}

void DefaultVisitor::visitTravelingWave(const TravelingWave &trw) {
    applyDefault(trw);
}

#ifdef ENABLE_OPAL_FEL
void DefaultVisitor::visitUndulator(const Undulator &u) {
    applyDefault(u);
}
#endif

void DefaultVisitor::visitVacuum(const Vacuum &vac) {
    applyDefault(vac);
}

void DefaultVisitor::visitVariableRFCavity(const VariableRFCavity &vcav) {
    applyDefault(vcav);
}

void DefaultVisitor::visitVariableRFCavityFringeField
                                (const VariableRFCavityFringeField &vcav) {
    applyDefault(vcav);
}

void DefaultVisitor::visitVerticalFFAMagnet(const VerticalFFAMagnet &mag) {
    applyDefault(mag);
}


void DefaultVisitor::visitBeamline(const Beamline &bl) {
    // Default behaviour: Apply algorithm to all beamline members.
    // If flip_local is true, track from right to left.
    bl.iterate(*this, local_flip);
}


void DefaultVisitor::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    if(fep.getReflectionFlag()) {
        local_flip = ! local_flip;
        fep.getElement()->accept(*this);
        local_flip = ! local_flip;
    } else {
        fep.getElement()->accept(*this);
    }
}


void DefaultVisitor::applyDefault(const ElementBase &)
{}
