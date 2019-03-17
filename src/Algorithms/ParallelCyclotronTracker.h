#ifndef OPAL_ParallelCyclotronTracker_HH
#define OPAL_ParallelCyclotronTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelCyclotronTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelCyclotron
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "Structure/DataSink.h"
#include "AbsBeamline/ElementBase.h"
#include <vector>
#include <tuple>

#include "Steppers/Steppers.h"

#include "Structure/MultiBunchDump.h"

template <class T, unsigned Dim>
class PartBunchBase;

// Class ParallelCyclotronTracker
// ------------------------------------------------------------------------

struct CavityCrossData {
    RFCavity * cavity;
    double sinAzimuth;
    double cosAzimuth;
    double perpenDistance;
};

class ParallelCyclotronTracker: public Tracker {

public:

    enum class MODE {
        UNDEFINED = -1,
        SINGLE = 0,
        SEO = 1,
        BUNCH = 2
    };

    // multi-bunch modes
    enum class MB_MODE {
        NONE   = 0,
        FORCE  = 1,
        AUTO   = 2
    };

    // multi-bunch binning type
    enum class MB_BINNING {
        GAMMA = 0,
        BUNCH = 1
    };

    typedef std::vector<double> dvector_t;
    typedef std::vector<int> ivector_t;
    typedef std::pair<double[8], Component *>      element_pair;
    typedef std::pair<ElementBase::ElementType, element_pair>        type_pair;
    typedef std::list<type_pair *>                 beamline_list;
    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    ParallelCyclotronTracker(const Beamline &bl, const PartData &data,
                             bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    ParallelCyclotronTracker(const Beamline &bl, PartBunchBase<double, 3> *bunch, DataSink &ds,
                             const PartData &data, bool revBeam,
                             bool revTrack, int maxSTEPS,
                             int timeIntegrator, int numBunch);

    virtual ~ParallelCyclotronTracker();

    /// Apply the algorithm to an Ring
    virtual void visitRing(const Ring &ring);

    /// Apply the algorithm to a Cyclotorn
    virtual void visitCyclotron(const Cyclotron &cycl);

    /// Apply the algorithm to a RFCavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to a Corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a Degrader
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a Diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a Drift.
    virtual void visitDrift(const Drift &);

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

    /// Apply the algorithm to a MultipoleT
    virtual void visitMultipoleT (const MultipoleT &);

    /// Apply the algorithm to a MultipoleTStraight
    virtual void visitMultipoleTStraight (const MultipoleTStraight &);

    /// Apply the algorithm to a MultipoleTCurvedConstRadius
    virtual void visitMultipoleTCurvedConstRadius (const MultipoleTCurvedConstRadius &);

    /// Apply the algorithm to a MultipoleTCurvedVarRadius
    virtual void visitMultipoleTCurvedVarRadius (const MultipoleTCurvedVarRadius &);

    /// Apply the algorithm to a Offset.
    virtual void visitOffset(const Offset &);

    /// Apply the algorithm to a Probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a RFQuadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply the algorithm to a SBend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a SBend3D.
    virtual void visitSBend3D(const SBend3D &);

    /// Apply the algorithm to a ScalingFFAGMagnet.
    virtual void visitScalingFFAGMagnet(const ScalingFFAGMagnet &bend);

    /// Apply the algorithm to a Separator.
    virtual void visitSeparator(const Separator &);

    /// Apply the algorithm to a Septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a Solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a charge stripper.
    virtual void visitStripper(const Stripper &);

    /// Apply the algorithm to a ParallelPlate, it is empty for cyclotrontracker .
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a CyclotronValley.it is empty for cyclotrontracker .
    virtual void visitCyclotronValley(const CyclotronValley &);

    /// Apply the algorithm to a VariabelRFCavity.
    virtual void visitVariableRFCavity(const VariableRFCavity &cav);

    /// Apply the algorithm to a VariabelRFCavity.
    virtual void visitVariableRFCavityFringeField
                                      (const VariableRFCavityFringeField &cav);

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// set total number of tracked bunches
    inline void setNumBunch(int n) { numBunch_m = n; }

    /// get total number of tracked bunches
    inline int  getNumBunch() { return numBunch_m; }

    /// set the working sub-mode for multi-bunch mode: "FORCE" or "AUTO"
    void setMultiBunchMode(const std::string& mbmode);

    // set binning type
    void setMultiBunchBinning(std::string binning);

    /// set the scale for binning in multi-bunch mode
    void setMultiBunchEta(const double& eta) { eta_m = eta; };

    /// set last dumped step
    inline void setLastDumpedStep(const int para) {lastDumpedStep_m = para ; }

    /// set the control parameter for "AUTO" sub-mode
    inline void setParaAutoMode(const double para) {CoeffDBunches_m = para; }

    ///@{ Method for restart
    inline void setPr(double x) {referencePr = x;}
    inline void setPt(double x) {referencePt = x;}
    inline void setPz(double x) {referencePz = x;}
    inline void setR(double x) {referenceR = x;}
    inline void setTheta(double x) {referenceTheta = x;}
    inline void setZ(double x) {referenceZ = x;}
    inline void setBeGa(double x) {bega = x;}
    inline void setPhi(double x) {referencePhi = x;}
    inline void setPsi(double x) {referencePsi = x;}
    inline void setPreviousH5Local(bool x) {previousH5Local = x;}
    ///@}
    void bgf_main_collision_test();
    void initializeBoundaryGeometry();

private:

    // Not implemented.
    ParallelCyclotronTracker();
    ParallelCyclotronTracker(const ParallelCyclotronTracker &);
    void operator=(const ParallelCyclotronTracker &);

    beamline_list FieldDimensions;
    std::list<Component *> myElements;
    Beamline *itsBeamline;

    DataSink *itsDataSink;

    std::unique_ptr<MultiBunchDump> itsMBDump_m;

    BoundaryGeometry *bgf_m;

    /// The maximal number of steps the system is integrated
    int maxSteps_m;

    /// The positive axes unit vectors
    static Vector_t const xaxis;
    static Vector_t const yaxis;
    static Vector_t const zaxis;

    /// The scale factor for dimensionless variables
    double scaleFactor_m;

    /// The reference variables
    double bega;
    double referenceR;
    double referenceTheta;
    double referenceZ = 0.0;

    double referencePr;
    double referencePt;
    double referencePz = 0.0;
    double referencePtot;

    double referencePsi;
    double referencePhi;

    bool spiral_flag = false;

    Vector_t PreviousMeanP;

    bool previousH5Local;

    double sinRefTheta_m;
    double cosRefTheta_m;

    /// The number of bunches specified in TURNS of RUN commond
    int numBunch_m;

    // 0 for single bunch (default),
    // 1 for FORCE,
    // 2 for AUTO
    MB_MODE multiBunchMode_m;
    
    // 0 for GAMMA (default),
    // 1 for BUNCH
    MB_BINNING binningType_m;

    // control parameter for AUTO multi-bunch mode
    double CoeffDBunches_m;

    int lastDumpedStep_m;

    double PathLength_m;

    void MtsTracker();

    void GenericTracker();
    bool getFieldsAtPoint(const double &t, const size_t &Pindex, Vector_t &Efield, Vector_t &Bfield);

    /*
      Local Variables both used by the integration methods
    */

    long long step_m;
    long long restartStep0_m;

    int turnnumber_m;

    double eta_m; // parameter for reset bin in multi-bunch run

    // record how many bunches has already been injected. ONLY FOR MPM
    int BunchCount_m;

    // used for automatic injection in multi-bunch mode
    double RLastTurn_m , RThisTurn_m;

    // external field arrays for dumping
    Vector_t FDext_m[2], extE_m, extB_m;

    const int myNode_m;
    const size_t initialLocalNum_m;
    const size_t initialTotalNum_m;

    std::ofstream outfTheta0_m;
    std::ofstream outfTheta1_m;
    std::ofstream outfTheta2_m;
    std::ofstream outfThetaEachTurn_m;


    //store the data of the beam which are required for injecting a new bunch for multibunch
    /// filename
    std::string onebunch_m;

    void openFiles(std::string fn);
    void closeFiles();

    void buildupFieldList(double BcParameter[], ElementBase::ElementType elementType, Component *elptr);

    // angle range [0~2PI) degree
    double calculateAngle(double x, double y);
    // angle range [-PI~PI) degree
    double calculateAngle2(double x, double y);

    bool readOneBunchFromFile(const size_t BeamCount);
    void saveOneBunch();
    
    void updateParticleBins_m();

    bool checkGapCross(Vector_t Rold, Vector_t Rnew, RFCavity * rfcavity, double &DistOld);
    bool RFkick(RFCavity * rfcavity, const double t, const double dt, const int Pindex);

    bool getTunes(dvector_t &t,  dvector_t &r,  dvector_t &z, int lastTurn, double &nur, double &nuz);

    IpplTimings::TimerRef IntegrationTimer_m;
    IpplTimings::TimerRef DumpTimer_m ;
    IpplTimings::TimerRef TransformTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;

    Vector_t calcMeanR() const;

    Vector_t calcMeanP() const;

    void repartition(); // Do repartition between nodes if step_m is multiple of Options::repartFreq

    // Transform the x- and y-parts of a particle attribute (position, momentum, fields) from the
    // global reference frame to the local reference frame.

    // phi is the angle of the bunch measured counter-clockwise from the positive x-axis.
    void globalToLocal(ParticleAttrib<Vector_t> & vectorArray,
                       double phi, Vector_t const translationToGlobal = 0);

    // Transform the x- and y-parts of a particle attribute (position, momentum, fields) from the
    // local reference frame to the global reference frame.
    void localToGlobal(ParticleAttrib<Vector_t> & vectorArray,
                       double phi, Vector_t const translationToGlobal = 0);

    // Overloaded version of globalToLocal using a quaternion for 3D rotation
    inline void globalToLocal(ParticleAttrib<Vector_t> & vectorArray,
                              Quaternion_t const quaternion,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using a quaternion for 3D rotation
    inline void localToGlobal(ParticleAttrib<Vector_t> & vectorArray,
                              Quaternion_t const quaternion,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of globalToLocal using phi and theta for pseudo 3D rotation
    inline void globalToLocal(ParticleAttrib<Vector_t> & particleVectors,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using phi and theta for pseudo 3D rotation
    inline void localToGlobal(ParticleAttrib<Vector_t> & particleVectors,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of globalToLocal using phi and theta for pseudo 3D rotation, single vector
    inline void globalToLocal(Vector_t & myVector,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using phi and theta for pseudo 3D rotation, single vector
    inline void localToGlobal(Vector_t & myVector,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Rotate the particles by an angle and axis defined in the quaternion (4-vector w,x,y,z)
    inline void rotateWithQuaternion(ParticleAttrib<Vector_t> & vectorArray, Quaternion_t const quaternion);

    // Returns the quaternion of the rotation from vector u to vector v
    inline void getQuaternionTwoVectors(Vector_t u, Vector_t v, Quaternion_t & quaternion);

    // Normalization of a quaternion
    inline void normalizeQuaternion(Quaternion_t & quaternion);

    // Normalization of a quaternion
    inline void normalizeVector(Vector_t & vector);

    // Some rotations
    inline void rotateAroundZ(ParticleAttrib<Vector_t> & particleVectors, double const phi);
    inline void rotateAroundX(ParticleAttrib<Vector_t> & particleVectors, double const psi);
    inline void rotateAroundZ(Vector_t & myVector, double const phi);
    inline void rotateAroundX(Vector_t & myVector, double const psi);

    // Push particles for time h.
    // Apply effects of RF Gap Crossings.
    // Update time and path length.
    // Unit assumptions: [itsBunch_m->R] = m, [itsBunch_m->P] = 1, [h] = s, [c] = m/s, [itsBunch_m->getT()] = s
    void push(double h);

    // Kick particles for time h
    // The fields itsBunch_m->Bf, itsBunch_m->Ef are used to calculate the forces
    void kick(double h);

    // Apply the trilogy half push - kick - half push,
    // considering only external fields
    void borisExternalFields(double h);

    // apply the plugin elements: probe, collimator, stripper, septum
    void applyPluginElements(const double dt);

    // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global apeture
    bool deleteParticle();

    std::ofstream outfTrackOrbit_m;

    void initTrackOrbitFile();

    void singleParticleDump();

    void bunchDumpStatData();
    
    // @param azimuth (global) [deg] of dump
    void bunchDumpStatDataPerBin(const double& azimuth);

    void bunchDumpPhaseSpaceData();

    void evaluateSpaceChargeField();

    void initDistInGlobalFrame();

    void checkNumPart(std::string s);

    // we store a pointer explicitly to the Ring
    Ring* opalRing_m;


    // If Ring is defined take the harmonic number from Ring; else use
    // cyclotron
    double getHarmonicNumber() const;

    typedef std::function<bool(const double&,
                               const size_t&,
                               Vector_t&,
                               Vector_t&)> function_t;

    std::unique_ptr< Stepper<function_t> > itsStepper_mp;

    struct settings {
        int scSolveFreq;
        int stepsPerTurn;
        double azimuth_angle0;
        double azimuth_angle1;
        double azimuth_angle2;
        double deltaTheta;
        int stepsNextCheck;
    } setup_m;

    MODE mode_m;

    stepper::INTEGRATOR stepper_m;

    void update_m(double& t, const double& dt, const bool& dumpEachTurn);

    /*!
     * @returns the time t [ns], time step dt [ns] and the azimuth angle [rad]
     */
    std::tuple<double, double, double> initializeTracking_m();

    void finalizeTracking_m(dvector_t& Ttime,
                            dvector_t& Tdeltr,
                            dvector_t& Tdeltz,
                            ivector_t& TturnNumber);

    void seoMode_m(double& t, const double dt, bool& dumpEachTurn,
                   dvector_t& Ttime, dvector_t& Tdeltr,
                   dvector_t& Tdeltz, ivector_t& TturnNumber);

    void singleMode_m(double& t, const double dt, bool& dumpEachTurn, double& oldReferenceTheta);

    void bunchMode_m(double& t, const double dt, bool& dumpEachTurn);

    void gapCrossKick_m(size_t i, double t, double dt,
                        const Vector_t& Rold, const Vector_t& Pold);


    inline void dumpAzimuthAngles_m(const double& t,
                                    const Vector_t& R,
                                    const Vector_t& P,
                                    const double& oldReferenceTheta,
                                    const double& temp_meanTheta);

    inline void dumpThetaEachTurn_m(const double& t,
                                    const Vector_t& R,
                                    const Vector_t& P,
                                    const double& temp_meanTheta,
                                    bool& dumpEachTurn);

    void computeSpaceChargeFields_m();

    bool computeExternalFields_m(const size_t& i,
                                 const double& t,
                                 Vector_t& Efield,
                                 Vector_t& Bfield);

    void injectBunch_m(bool& flagTransition);

};

/**
 * @param x
 * @param y
 *
 * @return angle range [0~2PI) degree
 */
inline
double ParallelCyclotronTracker::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    if (thetaXY < 0) return thetaXY + Physics::two_pi;
    return thetaXY;
}

/**
 * @param x
 * @param y
 *
 * @return angle range [-PI~PI) degree
 */
inline
double ParallelCyclotronTracker::calculateAngle2(double x, double y)
{ return atan2(y,x); }

#endif // OPAL_ParallelCyclotronTracker_HH
