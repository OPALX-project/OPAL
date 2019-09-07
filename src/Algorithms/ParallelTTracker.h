#ifndef OPAL_ParallelTTracker_HH
#define OPAL_ParallelTTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelTTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelTTracker
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "Steppers/BorisPusher.h"
#include "Structure/DataSink.h"
#include "Algorithms/StepSizeConfig.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/OrbitThreader.h"
#include "Algorithms/IndexMap.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"
#include "Solvers/WakeFunction.hh"

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

    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);


    /// Apply the algorithm to a Corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a Degrader.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a Diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a Drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a Undulator.
    virtual void visitUndulator(const Undulator &);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a Lambertson.
    virtual void visitLambertson(const Lambertson &);

    /// Apply the algorithm to a Marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a Monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a Multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to a MultipoleT.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to a Probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RFCavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a RFCavity.
    virtual void visitTravelingWave(const TravelingWave &);

    /// Apply the algorithm to a RFQuadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply the algorithm to a SBend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a Separator.
    virtual void visitSeparator(const Separator &);

    /// Apply the algorithm to a Septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a Solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a Solenoid.
    virtual void visitSource(const Source &);

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);


    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

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

    unsigned long totalParticlesInSimulation_m;

#ifdef OPAL_DKS
    DKSOPAL *dksbase;
    void *r_ptr;
    void *p_ptr;
    void *Ef_ptr;
    void *Bf_ptr;
    void *dt_ptr;

    int stream1;
    int stream2;

    int numDeviceElements;

    void setupDKS();
    void allocateDeviceMemory();
    void freeDeviceMemory();
    void sendRPdt();
    void sendEfBf();
    void pushParticlesDKS(bool send = true);
    void kickParticlesDKS();

#endif

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
    void computeUndulator(IndexMap::value_t &elements);
    void computeSpaceChargeFields(unsigned long long step);
    // void prepareOpalBeamlineSections();
    void dumpStats(long long step, bool psDump, bool statDump);
    void setOptionalVariables();
    bool hasEndOfLineReached();
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

inline void ParallelTTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, itsBunch_m);
}

inline void ParallelTTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, itsBunch_m);
}


inline void ParallelTTracker::visitCCollimator(const CCollimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}


inline void ParallelTTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch_m);
}


inline void ParallelTTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, itsBunch_m);
}


inline void ParallelTTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, itsBunch_m);
}


inline void ParallelTTracker::visitUndulator(const Undulator &u) {
    itsOpalBeamline_m.visit(u, *this, itsBunch_m);
}

inline void ParallelTTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}


inline void ParallelTTracker::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}


inline void ParallelTTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, itsBunch_m);
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

inline void ParallelTTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}


inline void ParallelTTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, itsBunch_m);
}

inline void ParallelTTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}


inline void ParallelTTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, itsBunch_m);
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

inline void ParallelTTracker::visitParallelPlate(const ParallelPlate &pplate) {
    itsOpalBeamline_m.visit(pplate, *this, itsBunch_m);
}

inline void ParallelTTracker::visitCyclotronValley(const CyclotronValley &cv) {
    itsOpalBeamline_m.visit(cv, *this, itsBunch_m);
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

#ifdef OPAL_DKS
inline void ParallelTTracker::setupDKS() {
    dksbase = new DKSOPAL;
    dksbase->setAPI("Cuda");
    dksbase->setDevice("-gpu");
    dksbase->initDevice();

    dksbase->createStream(stream1);
    dksbase->createStream(stream2);
}


inline void ParallelTTracker::allocateDeviceMemory() {

    int ierr;
    r_ptr = dksbase->allocateMemory<Vector_t>(itsBunch_m->getLocalNum(), ierr);
    p_ptr = dksbase->allocateMemory<Vector_t>(itsBunch_m->getLocalNum(), ierr);
    Ef_ptr = dksbase->allocateMemory<Vector_t>(itsBunch_m->getLocalNum(), ierr);
    Bf_ptr = dksbase->allocateMemory<Vector_t>(itsBunch_m->getLocalNum(), ierr);
    dt_ptr = dksbase->allocateMemory<double>(itsBunch_m->getLocalNum(), ierr);

    if (Ippl::getNodes() == 1) {
        dksbase->registerHostMemory(&itsBunch_m->R[0], itsBunch_m->getLocalNum());
        dksbase->registerHostMemory(&itsBunch_m->P[0], itsBunch_m->getLocalNum());
        dksbase->registerHostMemory(&itsBunch_m->Ef[0], itsBunch_m->getLocalNum());
        dksbase->registerHostMemory(&itsBunch_m->Bf[0], itsBunch_m->getLocalNum());
        dksbase->registerHostMemory(&itsBunch_m->dt[0], itsBunch_m->getLocalNum());

    }

    numDeviceElements = itsBunch_m->getLocalNum();
}

inline  void ParallelTTracker::freeDeviceMemory() {
    dksbase->freeMemory<Vector_t>(r_ptr, numDeviceElements);
    dksbase->freeMemory<Vector_t>(p_ptr, numDeviceElements);
    dksbase->freeMemory<Vector_t>(Ef_ptr, numDeviceElements);
    dksbase->freeMemory<Vector_t>(Bf_ptr, numDeviceElements);
    dksbase->freeMemory<double>(dt_ptr, numDeviceElements);

    if (Ippl::getNodes() == 1) {
        dksbase->unregisterHostMemory(&itsBunch_m->R[0]);
        dksbase->unregisterHostMemory(&itsBunch_m->P[0]);
        dksbase->unregisterHostMemory(&itsBunch_m->Ef[0]);
        dksbase->unregisterHostMemory(&itsBunch_m->Bf[0]);
        dksbase->unregisterHostMemory(&itsBunch_m->dt[0]);
    }
}

inline void ParallelTTracker::sendRPdt() {
    dksbase->writeDataAsync<Vector_t>(r_ptr, &itsBunch_m->R[0], itsBunch_m->getLocalNum(), stream1);
    dksbase->writeDataAsync<Vector_t>(p_ptr, &itsBunch_m->P[0], itsBunch_m->getLocalNum(), stream1);
    dksbase->writeDataAsync<double>(dt_ptr, &itsBunch_m->dt[0], itsBunch_m->getLocalNum(), stream1);
}

inline void ParallelTTracker::sendEfBf() {
    dksbase->writeDataAsync<Vector_t>(Ef_ptr, &itsBunch_m->Ef[0],
                                      itsBunch_m->getLocalNum(), stream1);
    dksbase->writeDataAsync<Vector_t>(Bf_ptr, &itsBunch_m->Bf[0],
                                      itsBunch_m->getLocalNum(), stream1);
}

inline void ParallelTTracker::pushParticlesDKS(bool send) {

    //send data to the GPU
    if (send)
        sendRPdt();
    //execute particle push on the GPU
    dksbase->callParallelTTrackerPush(r_ptr, p_ptr, dt_ptr, itsBunch_m->getLocalNum(),
                                      Physics::c, stream1);
    //get particles back to CPU
    dksbase->readData<Vector_t>(r_ptr, &itsBunch_m->R[0], itsBunch_m->getLocalNum(), stream1);
}

inline
void ParallelTTracker::kickParticlesDKS() {
    //send data to the GPU
    sendEfBf();
    sendRPdt();

    //sync device
    dksbase->syncDevice();

    //execute kick on the GPU
    dksbase->callParallelTTrackerKick(r_ptr, p_ptr, Ef_ptr, Bf_ptr, dt_ptr,
                                      itsReference.getQ(), itsReference.getM(),
                                      itsBunch_m->getLocalNum(), Physics::c, stream2);

    dksbase->syncDevice();

    //get data back from GPU
    dksbase->readDataAsync<Vector_t>(p_ptr, &itsBunch_m->P[0], itsBunch_m->getLocalNum(), stream2);

}
#endif

#endif // OPAL_ParallelTTracker_HH