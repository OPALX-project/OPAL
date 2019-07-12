#ifndef OPAL_ParallelSliceTracker_HH
#define OPAL_ParallelSliceTracker_HH

#include <list>
#include <memory>

#include "Algorithms/bet/EnvelopeBunch.h"

#include "Algorithms/Vektor.h"
#include "Algorithms/Tracker.h"
#include "Structure/DataSink.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

#include "Algorithms/OrbitThreader.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"

#include "Algorithms/ParallelTTracker.h"

#include "Beamlines/Beamline.h"
#include "Elements/OpalBeamline.h"

class ParallelSliceTracker: public Tracker {

public:

    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelSliceTracker(const Beamline &bl, const PartData &data,
                                  bool revBeam, bool revTrack);

    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ParallelSliceTracker(const Beamline &bl,
                                  EnvelopeBunch &bunch,
                                  DataSink &ds,
                                  const PartData &data,
                                  bool revBeam,
                                  bool revTrack,
                                  const std::vector<unsigned long long> &maxSTEPS,
                                  double zstart,
                                  const std::vector<double> &zstop,
                                  const std::vector<double> &dt);

    virtual ~ParallelSliceTracker();

    virtual void visitAlignWrapper(const AlignWrapper &);
    virtual void visitBeamBeam(const BeamBeam &);
    virtual void visitCCollimator(const CCollimator &);
    virtual void visitCorrector(const Corrector &);
    virtual void visitDegrader(const Degrader &);
    virtual void visitDiagnostic(const Diagnostic &);
    virtual void visitDrift(const Drift &);
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);
    virtual void visitLambertson(const Lambertson &);
    virtual void visitMarker(const Marker &);
    virtual void visitMonitor(const Monitor &);
    virtual void visitMultipole(const Multipole &);
    virtual void visitProbe(const Probe &);
    virtual void visitRBend(const RBend &);
    virtual void visitRFCavity(const RFCavity &);
    virtual void visitTravelingWave(const TravelingWave &);
    virtual void visitRFQuadrupole(const RFQuadrupole &);
    virtual void visitSBend(const SBend &);
    virtual void visitSeparator(const Separator &);
    virtual void visitSeptum(const Septum &);
    virtual void visitSolenoid(const Solenoid &);
    virtual void visitSource(const Source &);
    virtual void visitParallelPlate(const ParallelPlate &);
    virtual void visitCyclotronValley(const CyclotronValley &);

    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);


private:

    ParallelSliceTracker();
    ParallelSliceTracker(const ParallelSliceTracker &);
    void operator=(const ParallelSliceTracker &);

    void saveCavityPhases();
    void restoreCavityPhases();
    void updateRFElement(std::string elName, double maxPhi);
    void printRFPhases();

    void findStartPosition(const BorisPusher &pusher);
    void autophaseCavities(const BorisPusher &pusher);

    double currentSimulationTime_m;
    bool globalEOL_m;

    OpalBeamline itsOpalBeamline_m;

    EnvelopeBunch *itsBunch_m;

    DataSink *itsDataSink_m;

    Vector_t RefPartR_m;
    Vector_t RefPartP_m;

    double pathLength_m;

    /// where to start
    double zstart_m;

    /// where to stop
    std::queue<double> zStop_m;

    double dtCurrentTrack_m;
    std::queue<double> dtAllTracks_m;

    /// The maximal number of steps the system is integrated per TRACK
    std::queue<unsigned long long> localTrackSteps_m;

    CoordinateSystemTrafo referenceToLabCSTrafo_m;

    // Fringe fields for entrance and exit of magnetic elements.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

    void FieldsOutput(double z, double Ex, double Ey, double Ez,
                      double Bx, double By, double Bz);

    void kickParticles();

    void updateReferenceParticle(const BorisPusher &pusher);
    // void updateReferenceParticle();
    void updateSpaceOrientation(const bool &move = false);
    void kickReferenceParticle(const Vector_t &externalE,
                               const Vector_t &externalB);
    void writePhaseSpace(const long long step, const double &sposRef);
    void writeLastStepPhaseSpace(const long long step, const double &sposRef);

    void prepareSections();
    void handleAutoPhasing();
    void timeIntegration();
    void computeExternalFields(OrbitThreader &oth);
    void switchElements(double scaleMargin = 3.0);
    void computeSpaceChargeFields();
    void dumpStats(long long step);
    bool hasEndOfLineReached();
    void handleRestartRun();
    void setTime();
    void setLastStep();
    unsigned long long getMaxSteps(std::queue<unsigned long long> numSteps);
    void selectDT();
    void changeDT();
    IpplTimings::TimerRef timeIntegrationTimer1_m;
    IpplTimings::TimerRef timeIntegrationTimer2_m;
    IpplTimings::TimerRef timeFieldEvaluation_m ;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef WakeFieldTimer_m;
};

inline void ParallelSliceTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    itsOpalBeamline_m.visit(wrap, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitBeamBeam(const BeamBeam &bb) {
    itsOpalBeamline_m.visit(bb, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitCCollimator(const CCollimator &coll) {
    // itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    // itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitCorrector(const Corrector &corr) {
    itsOpalBeamline_m.visit(corr, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitDegrader(const Degrader &deg) {
    itsOpalBeamline_m.visit(deg, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitDiagnostic(const Diagnostic &diag) {
    itsOpalBeamline_m.visit(diag, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitLambertson(const Lambertson &lamb) {
    itsOpalBeamline_m.visit(lamb, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitMarker(const Marker &marker) {
    itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitMonitor(const Monitor &mon) {
    itsOpalBeamline_m.visit(mon, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitProbe(const Probe &prob) {
    itsOpalBeamline_m.visit(prob, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitRBend(const RBend &bend) {
    // itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitRFCavity(const RFCavity &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitTravelingWave(const TravelingWave &as) {
    itsOpalBeamline_m.visit(as, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    itsOpalBeamline_m.visit(rfq, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitSBend(const SBend &bend) {
    // itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitSeparator(const Separator &sep) {
    itsOpalBeamline_m.visit(sep, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitSeptum(const Septum &sept) {
    itsOpalBeamline_m.visit(sept, *this, itsBunch_m);
}


inline void ParallelSliceTracker::visitSolenoid(const Solenoid &solenoid) {
    itsOpalBeamline_m.visit(solenoid, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitSource(const Source &source) {
    itsOpalBeamline_m.visit(source, *this, itsBunch_m);
}

inline void ParallelSliceTracker::visitParallelPlate(const ParallelPlate &pplate) {
    //do nothing.
}

inline void ParallelSliceTracker::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}

inline void ParallelSliceTracker::kickParticles() {
}

inline void ParallelSliceTracker::updateSpaceOrientation(const bool &move) {
    itsBunch_m->calcBeamParameters();
}

inline void ParallelSliceTracker::kickReferenceParticle(const Vector_t &externalE, const Vector_t &externalB) {
}

inline void ParallelSliceTracker::writePhaseSpace(const long long step, const double &sposRef) {
    Inform msg("ParallelSliceTracker");
    Vector_t externalE, externalB;
    Vector_t FDext[6];  // = {BHead, EHead, BRef, ERef, BTail, ETail}.

    /*
     * Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the
     * centroid location. We are sampling the electric and magnetic fields at
     * the back, front and center of the beam.
     */
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);

    Vector_t pos[3];
    pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
    pos[1] = Vector_t(0.0, 0.0, sposRef);
    pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

    for(int k = 0; k < 3; ++k) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        double bunch_time = itsBunch_m->getT() - 0.5 * itsBunch_m->getdT();
        itsOpalBeamline_m.getFieldAt(pos[k], itsBunch_m->get_rmean(),
                                      bunch_time, externalE, externalB);

        FDext[2*k]   = externalB;
        FDext[2*k+1] = externalE * 1e-6;
    }

    if(step % Options::psDumpFreq == 0) {
        //itsDataSink->stashPhaseSpaceEnvelope(*itsBunch_m, FDext, rmax(2), sposRef, rmin(2));
        itsDataSink_m->dumpH5(*itsBunch_m, FDext, rmax(2), sposRef, rmin(2));
        msg     << *itsBunch_m << endl;
    }

    if(step % Options::statDumpFreq == 0) {
        itsDataSink_m->dumpSDDS(*itsBunch_m, FDext, rmax(2), sposRef, rmin(2));
    }
}

inline void ParallelSliceTracker::writeLastStepPhaseSpace(const long long step, const double &sposRef) {
    Inform msg("ParallelSliceTracker");
    if(itsBunch_m->isValid_m) {
        Vector_t externalE, externalB;
        Vector_t FDext[6];
        Vector_t rmin, rmax;
        itsBunch_m->get_bounds(rmin, rmax);

        Vector_t pos[3];
        pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
        pos[1] = Vector_t(0.0, 0.0, sposRef);
        pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

        for(int k = 0; k < 3; ++k) {
            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            FDext[2*k]   = externalB;
            FDext[2*k+1] = externalE * 1e-6;
        }

        itsDataSink_m->dumpSDDS(*itsBunch_m, FDext, rmax(2), sposRef, rmin(2));
    } else
        INFOMSG("* Invalid bunch! No statistics dumped." << endl);
}

inline unsigned long long ParallelSliceTracker::getMaxSteps(std::queue<unsigned long long> numSteps) {
    unsigned long long totalNumSteps = 0;

    while (numSteps.size() > 0) {
        totalNumSteps += numSteps.front();
        numSteps.pop();
    }

    return totalNumSteps;
}


#endif