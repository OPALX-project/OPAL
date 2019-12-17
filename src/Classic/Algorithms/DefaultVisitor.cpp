// ------------------------------------------------------------------------
// $RCSfile: DefaultVisitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DefaultVisitor
//   Defines the default interface for a BeamlineVisitor.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:04 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/BeamStripping.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/Patch.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/VerticalFFAMagnet.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/Source.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/Stripper.h"

#include "Algorithms/MapIntegrator.h"
#include "Algorithms/TrackIntegrator.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

#include "ComponentWrappers/CorrectorWrapper.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "ComponentWrappers/SBendWrapper.h"
#include "ComponentWrappers/CyclotronWrapper.h"

#include "AbsBeamline/Ring.h" // OPAL file

// Class DefaultVisitor
// ------------------------------------------------------------------------

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


void DefaultVisitor::visitBeamBeam(const BeamBeam &bb) {
    applyDefault(bb);
}

void DefaultVisitor::visitBeamStripping(const BeamStripping &bstp) {
    applyDefault(bstp);
}

void DefaultVisitor::visitCCollimator(const CCollimator &coll) {
    applyDefault(coll);
}

void DefaultVisitor::visitCyclotron(const Cyclotron &cyc) {
    applyDefault(cyc);
}

void DefaultVisitor::visitComponent(const Component &comp) {
    applyDefault(comp);
}

void DefaultVisitor::visitCorrector(const Corrector &corr) {
    applyDefault(corr);
}

void DefaultVisitor::visitDegrader(const Degrader &deg) {
    applyDefault(deg);
}

void DefaultVisitor::visitDiagnostic(const Diagnostic &diag) {
    applyDefault(diag);
}

void DefaultVisitor::visitDrift(const Drift &drf) {
    applyDefault(drf);
}

void DefaultVisitor::visitUndulator(const Undulator &u) {
    applyDefault(u);
}


void DefaultVisitor::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDefault(coll);
}

void DefaultVisitor::visitLambertson(const Lambertson &lamb) {
    applyDefault(lamb);
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

void DefaultVisitor::visitRing(const Ring &ring) {
   applyDefault(ring);
}


void DefaultVisitor::visitPatch(const Patch &patch) {
    applyDefault(patch);
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

void DefaultVisitor::visitVariableRFCavity(const VariableRFCavity &vcav) {
    applyDefault(vcav);
}

void DefaultVisitor::visitVariableRFCavityFringeField
                                (const VariableRFCavityFringeField &vcav) {
    applyDefault(vcav);
}

void DefaultVisitor::visitRFCavity(const RFCavity &cav) {
    applyDefault(cav);
}


void DefaultVisitor::visitTravelingWave(const TravelingWave &trw) {
    applyDefault(trw);
}


void DefaultVisitor::visitRFQuadrupole(const RFQuadrupole &quad) {
    applyDefault(quad);
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

void DefaultVisitor::visitVerticalFFAMagnet(const VerticalFFAMagnet &mag) {
    applyDefault(mag);
}

void DefaultVisitor::visitSeparator(const Separator &sep) {
    applyDefault(sep);
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


void DefaultVisitor::visitParallelPlate(const ParallelPlate &pplate) {
    applyDefault(pplate);
}

void DefaultVisitor::visitCyclotronValley(const CyclotronValley &cv) {
    applyDefault(cv);
}

void DefaultVisitor::visitStripper(const Stripper &cv) {
    applyDefault(cv);
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


void DefaultVisitor::visitAlignWrapper(const AlignWrapper &wrap) {
    // Default behaviour: Apply algorithm to the non-offset element.
    wrap.getElement()->accept(*this);
}


void DefaultVisitor::visitCorrectorWrapper(const CorrectorWrapper &wrap) {
    visitCorrector(wrap);
}

void DefaultVisitor::visitCyclotronWrapper(const CyclotronWrapper &wrap) {
    visitCyclotron(wrap);
}


void DefaultVisitor::visitMultipoleWrapper(const MultipoleWrapper &wrap) {
    visitMultipole(wrap);
}


void DefaultVisitor::visitRBendWrapper(const RBendWrapper &wrap) {
    visitRBend(wrap);
}


void DefaultVisitor::visitSBendWrapper(const SBendWrapper &wrap) {
    visitSBend(wrap);
}


void DefaultVisitor::visitIntegrator(const Integrator &i) {
    // Default: cannot use integrator.
    i.getElement()->accept(*this);
}


void DefaultVisitor::visitMapIntegrator(const MapIntegrator &i) {
    // Default: cannot use integrator.
    i.getElement()->accept(*this);
}


void DefaultVisitor::visitTrackIntegrator(const TrackIntegrator &i) {
    // Default: cannot use integrator.
    i.getElement()->accept(*this);
}


void DefaultVisitor::applyDefault(const ElementBase &)
{}
