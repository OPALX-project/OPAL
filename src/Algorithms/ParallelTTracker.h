//
// Class ParallelTTracker
//   OPAL-T tracker.
//   The visitor class for tracking particles with time as independent
//   variable.
//
// Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020, Christof Metzger-Kraus
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
#ifndef OPAL_ParallelTTracker_HH
#define OPAL_ParallelTTracker_HH

#include "Algorithms/Tracker.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"
#include "Algorithms/StepSizeConfig.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/OrbitThreader.h"
#include "Algorithms/IndexMap.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"
#ifdef ENABLE_OPAL_FEL
#include "AbsBeamline/Undulator.h"
#endif
#include "AbsBeamline/Vacuum.h"

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Solvers/WakeFunction.h"

#include <list>
#include <vector>

class ParticleMatterInteractionHandler;

class ParallelTTracker: public Tracker {

public:
    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl,
                              const PartData &data,
                              bool revBeam,
                              bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelTTracker(const Beamline &bl,
                              PartBunchBase<double, 3> *bunch,
                              DataSink &ds,
                              const PartData &data,
                              bool revBeam,
                              bool revTrack,
                              const std::vector<unsigned long long> &maxSTEPS,
                              double zstart,
                              const std::vector<double> &zstop,
                              const std::vector<double> &dt);


    virtual ~ParallelTTracker();

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to a closed orbit corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a degrader.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a source.
    virtual void visitSource(const Source &);

    /// Apply the algorithm to a traveling wave.
    virtual void visitTravelingWave(const TravelingWave &);

#ifdef ENABLE_OPAL_FEL
    /// Apply the algorithm to an undulator.
    virtual void visitUndulator(const Undulator &);
#endif

    /// Apply the algorithm to a vacuum space.
    virtual void visitVacuum(const Vacuum &);

private:

    // Not implemented.
    ParallelTTracker();
    ParallelTTracker(const ParallelTTracker &);
    void operator=(const ParallelTTracker &);

    /******************** STATE VARIABLES ***********************************/

    DataSink *itsDataSink_m;

    OpalBeamline itsOpalBeamline_m;

    bool globalEOL_m;

    bool wakeStatus_m;

    bool deletedParticles_m;

    WakeFunction* wakeFunction_m;

    double pathLength_m;

    /// where to start
    double zstart_m;

    /// stores informations where to change the time step and
    /// where to stop the simulation,
    /// the time step sizes and
    /// the number of time steps with each configuration
    StepSizeConfig stepSizes_m;

    double dtCurrentTrack_m;

    // This variable controls the minimal number of steps of emission (using bins)
    // before we can merge the bins
    int minStepforReBin_m;

    // The space charge solver crashes if we use less than ~10 particles.
    // This variable controls the number of particles to be emitted before we use
    // the space charge solver.
    size_t minBinEmitted_m;

    // this variable controls the minimal number of steps until we repartition the particles
    unsigned int repartFreq_m;

    unsigned int emissionSteps_m;

    size_t numParticlesInSimulation_m;

    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef fieldEvaluationTimer_m ;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;

    std::set<ParticleMatterInteractionHandler*> activeParticleMatterInteractionHandlers_m;
    bool particleMatterStatus_m;

    /********************** END VARIABLES ***********************************/

    void kickParticles(const BorisPusher &pusher);
    void pushParticles(const BorisPusher &pusher);
    void updateReferenceParticle(const BorisPusher &pusher);

    void writePhaseSpace(const long long step, bool psDump, bool statDump);

    /********** BEGIN AUTOPHSING STUFF **********/
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();
    void saveCavityPhases();
    void restoreCavityPhases();
    /************ END AUTOPHSING STUFF **********/

    void prepareSections();

    void timeIntegration1(BorisPusher & pusher);
    void timeIntegration2(BorisPusher & pusher);
    void selectDT(bool backTrack = false);
    void changeDT(bool backTrack = false);
    void emitParticles(long long step);
    void computeExternalFields(OrbitThreader &oth);
    void computeWakefield(IndexMap::value_t &elements);
    void computeParticleMatterInteraction(IndexMap::value_t elements, OrbitThreader &oth);
#ifdef ENABLE_OPAL_FEL
    void computeUndulator(IndexMap::value_t &elements);
#endif
    void computeSpaceChargeFields(unsigned long long step);
    // void prepareOpalBeamlineSections();
    void dumpStats(long long step, bool psDump, bool statDump);
    void setOptionalVariables();
    bool hasEndOfLineReached(const BoundingBox& globalBoundingBox);
    void handleRestartRun();
    void prepareEmission();
    void setTime();
    void doBinaryRepartition();

    void transformBunch(const CoordinateSystemTrafo &trafo);

    void updateReference(const BorisPusher &pusher);
    void updateRefToLabCSTrafo();
    void applyFractionalStep(const BorisPusher &pusher, double tau);
    void findStartPosition(const BorisPusher &pusher);
    void autophaseCavities(const BorisPusher &pusher);

    void evenlyDistributeParticles();
};


inline void ParallelTTracker::visitCCollimator(const CCollimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}

inline void ParallelTTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch_m);
}

inline void ParallelTTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, itsBunch_m);
}

inline void ParallelTTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}

inline void ParallelTTracker::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ParallelTTracker::visitMultipoleT(const MultipoleT &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}

inline void ParallelTTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, itsBunch_m);
}

inline void ParallelTTracker::visitRBend(const RBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}

inline void ParallelTTracker::visitRBend3D(const RBend3D &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}

inline void ParallelTTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ParallelTTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}

inline void ParallelTTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, itsBunch_m);
}

inline void ParallelTTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, itsBunch_m);
}

inline void ParallelTTracker::visitSource(const Source &source) {
    itsOpalBeamline_m.visit(source, *this, itsBunch_m);
}

inline void ParallelTTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

#ifdef ENABLE_OPAL_FEL
inline void ParallelTTracker::visitUndulator(const Undulator &u) {
    itsOpalBeamline_m.visit(u, *this, itsBunch_m);
}
#endif

inline void ParallelTTracker::visitVacuum(const Vacuum &vac) {
    itsOpalBeamline_m.visit(vac, *this, itsBunch_m);
}

inline void ParallelTTracker::kickParticles(const BorisPusher &pusher) {
    int localNum = itsBunch_m->getLocalNum();
    for (int i = 0; i < localNum; ++i)
        pusher.kick(itsBunch_m->R[i], itsBunch_m->P[i], itsBunch_m->Ef[i], itsBunch_m->Bf[i], itsBunch_m->dt[i]);
}

inline void ParallelTTracker::pushParticles(const BorisPusher &pusher) {
    itsBunch_m->switchToUnitlessPositions(true);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        pusher.push(itsBunch_m->R[i], itsBunch_m->P[i], itsBunch_m->dt[i]);
    }
    itsBunch_m->switchOffUnitlessPositions(true);
}

#endif // OPAL_ParallelTTracker_HH