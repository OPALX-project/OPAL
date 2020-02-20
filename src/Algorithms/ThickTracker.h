#ifndef OPAL_ThickTracker_HH
#define OPAL_ThickTracker_HH

// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickTracker
//
// ------------------------------------------------------------------------
//
// $Author: ganz_p $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "Structure/DataSink.h"

#include "Hamiltonian.h"

#include "MapAnalyser.h"

#include "Algorithms/IndexMap.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/BeamStripping.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RBend3D.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/TravelingWave.h"

#include "Elements/OpalBeamline.h"

#include "Structure/Beam.h"

//#include <array>
#include <cmath>

#include <tuple>

class BMultipoleField;

template <class T, unsigned Dim>
class PartBunchBase;



// Class ThickTracker
// ------------------------------------------------------------------------
/// Track using thick-lens algorithm.
// [p]
// Phase space coordinates numbering:
// [tab 3 b]
// [row]number [&]name          [&]unit  [/row]
// [row]0      [&]$x$           [&]metres [/row]
// [row]1      [&]$p_x/p_r$     [&]1      [/row]
// [row]2      [&]$y$           [&]metres [/row]
// [row]3      [&]$p_y/p_r$     [&]1      [/row]
// [row]4      [&]$v*delta_t$   [&]metres [/row]
// [row]5      [&]$delta_p/p_r$ [&]1      [/row]
// [/tab][p]
// Where $p_r$ is the constant reference momentum defining the reference
// frame velocity, $m$ is the rest mass of the particles, and $v$ is the
// instantaneous velocity of the particle.
// [p]
// Other units used:
// [tab 2 b]
// [row]quantity             [&]unit           [/row]
// [row]reference momentum   [&]electron-volts [/row]
// [row]velocity             [&]metres/second  [/row]
// [row]accelerating voltage [&]volts          [/row]
// [row]separator voltage    [&]volts          [/row]
// [row]frequencies          [&]hertz          [/row]
// [row]phase lags           [&]$2*pi$         [/row]
// [/tab][p]
// Approximations used:
// [ul]
// [li]
// [li]
// [li]
// [/ul]
//
// On going through an element, we use the following steps:
// To complete the map, we propagate the closed orbit and add that to the map.

class ThickTracker: public Tracker {

public:
    typedef Hamiltonian::series_t                       series_t;
    typedef MapAnalyser::fMatrix_t                      fMatrix_t;
    typedef MapAnalyser::cfMatrix_t                     cfMatrix_t;
    typedef FVps<double, 6>                             map_t;
    typedef FVector<double, 6>                          particle_t;
    typedef std::tuple<series_t, std::size_t, double>   tuple_t;
    typedef std::list<tuple_t>                          beamline_t;

public:

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ThickTracker(const Beamline &bl, const PartData &data,
                          bool revBeam, bool revTrack);

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ThickTracker(const Beamline &bl,
                          PartBunchBase<double, 3> *bunch,
                          Beam &beam,
                          DataSink &ds,
                          const PartData &data,
                          bool revBeam, bool revTrack,
                          const std::vector<unsigned long long> &maxSTEPS,
                          double zstart,
                          const std::vector<double> &zstop,
                          const std::vector<double> &dt,
                          const int& truncOrder);

    virtual ~ThickTracker();


    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a BeamStripping.
    virtual void visitBeamStripping(const BeamStripping &);

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
    /// Apply the algorithm to a Probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend(const RBend &);

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

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline &);

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    void prepareSections();

    /*
    void insertFringeField(SBend* pSBend, std::list<structMapTracking>& mBL, double& beta0,
             double& gamma0, double& P0, double& q, std::array<double,2>& entrFringe, std::string e);
     */
private:


    // Not implemented.
    ThickTracker() = delete;
    ThickTracker(const ThickTracker &)   = delete;
    void operator=(const ThickTracker &) = delete;

    void throwElementError_m(std::string element) {
        throw LogicalError("ThickTracker::execute()",
                           "Element '" + element + "' not supported.");
    }



    /*!
     * Tests the order of the elements in the beam line according to their position
     */
    void checkElementOrder_m();

    /*!
     * Inserts Drift maps in undefined beam line sections
     */
    void fillGaps_m();

    /*!
     * Tracks itsBunch_m trough beam line
     */
    void track_m();

    particle_t particleToVector_m(const Vector_t& R,
                                  const Vector_t& P) const;

    void vectorToParticle_m(const particle_t& particle,
                            Vector_t& R,
                            Vector_t& P) const;

    /*!
     * Advances itsBunch_m trough map
     * @param map Map of slice
     */
    void advanceParticles_m(const map_t& map);

    /*!
     * Applies map on particle
     * @param particle tracked particle
     * @param map Map of slice
     */
    void updateParticle_m(particle_t& particle,
                          const map_t& map);

    /*!
     * Dumps bunch in .stat or .h5 files
     */
    void dump_m();


    /*!
     * Updates itsBunch_m
     * @param spos position of tracking
     * @param step stepsize of applied map
     */
    void update_m(const double& spos,
                  const std::size_t& step);


    /*!
     * Writes map (Transfermap) in .map file
     * @map text for .map file
     */
    void write_m(const map_t& map);

    /*!
     * Concatenate map x and y
     * \f[
     *  y = x \circ y
     * \f]
     * @param x
     * @param y is result
     */
    void concatenateMaps_m(const map_t& x, map_t& y);


    /*!
     * Tracks Dispersion along beam line
     * writes it in .dispersion file.
     *
     *  \f[
     *  \begin{pmatrix} \eta_{x} \\ \eta_{p_x} \end{pmatrix}_{s_1}
     *  =
     *  \begin{pmatrix} R_{11} & R_{12} \\ R_{21} & R_{22} \end{pmatrix}
     *  \cdot
     *  \begin{pmatrix} \eta_{x} \\ \eta_{p_x} \end{pmatrix}_{s_0}
     *  +
     *  \begin{pmatrix} R_{16} \\ R_{26} \end{pmatrix}
     *  \f]
     *
     * @param tempMatrix accumulated Transfer map \f$R\f$ at pos
     * @param initialVal initial Dispersion { \f$\eta_{x0},\, \eta_{p_x0},\, \eta_{y0},\, \eta_{p_y0} \f$}
     * @param pos position of tracking
     */
    void advanceDispersion_m(fMatrix_t tempMatrix,
                             FMatrix<double, 1, 4> initialVal,
                             double pos);

    /*! :TODO:
     * Fringe fields for entrance of SBEND.
     *
     * @param edge
     * @param curve
     * @param field
     * @param scale
     */
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);

    /*! :TODO:
     * Fringe fields for exit of SBEND.
     *
     * @param edge
     * @param curve
     * @param field
     * @param scale
     */
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);


    Hamiltonian hamiltonian_m;
    MapAnalyser mapAnalyser_m;


    Vector_t RefPartR_m;
    Vector_t RefPartP_m;

    DataSink *itsDataSink_m;

    OpalBeamline itsOpalBeamline_m;

    double zstart_m;        ///< Start of beam line
    double zstop_m;         ///< End of beam line
    double threshold_m;     ///< Threshold for element overlaps and gaps
    beamline_t elements_m;  ///< elements in beam line

    CoordinateSystemTrafo referenceToLabCSTrafo_m;


    int truncOrder_m; ///< truncation order of map tracking

    IpplTimings::TimerRef mapCreation_m;    ///< creation of elements_m
    IpplTimings::TimerRef mapCombination_m; ///< map accumulation along elements_m -> Transfermap
    IpplTimings::TimerRef mapTracking_m;    ///< track particles trough maps of elements_m
};

inline void ThickTracker::visitAlignWrapper(const AlignWrapper &/*wrap*/) {
//     itsOpalBeamline_m.visit(wrap, *this, itsBunch_m);
    this->throwElementError_m("AlignWrapper");
}

inline void ThickTracker::visitBeamBeam(const BeamBeam &/*bb*/) {
//     itsOpalBeamline_m.visit(bb, *this, itsBunch_m);
    this->throwElementError_m("BeamBeam");
}

inline void ThickTracker::visitBeamStripping(const BeamStripping &bstp) {
    itsOpalBeamline_m.visit(bstp, *this, itsBunch_m);
}

inline void ThickTracker::visitCCollimator(const CCollimator &/*coll*/) {
//     itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
    this->throwElementError_m("CCollimator");
}


inline void ThickTracker::visitCorrector(const Corrector &/*coll*/) {
//     itsOpalBeamline_m.visit(corr, *this, itsBunch_m);
    this->throwElementError_m("Corrector");
}


inline void ThickTracker::visitDegrader(const Degrader &/*deg*/) {
//     itsOpalBeamline_m.visit(deg, *this, itsBunch_m);
    this->throwElementError_m("Degrader");
}


inline void ThickTracker::visitDiagnostic(const Diagnostic &/*diag*/) {
//     itsOpalBeamline_m.visit(diag, *this, itsBunch_m);
    this->throwElementError_m("Diagnostic");
}


inline void ThickTracker::visitDrift(const Drift &drift) {
    itsOpalBeamline_m.visit(drift, *this, itsBunch_m);


    double gamma = itsReference.getGamma();
    std::size_t nSlices = drift.getNSlices();
    double length       = drift.getElementLength();

    elements_m.push_back(std::make_tuple(hamiltonian_m.drift(gamma),
                                         nSlices,
                                         length));
}


inline void ThickTracker::visitFlexibleCollimator(const FlexibleCollimator &/*coll*/) {
//     itsOpalBeamline_m.visit(coll, *this, itsBunch_m);
    this->throwElementError_m("FlexibleCollimator");
}


inline void ThickTracker::visitLambertson(const Lambertson &/*lamb*/) {
//     itsOpalBeamline_m.visit(lamb, *this, itsBunch_m);
    this->throwElementError_m("Lambertson");
}


inline void ThickTracker::visitMarker(const Marker &/*marker*/) {
//     itsOpalBeamline_m.visit(marker, *this, itsBunch_m);
//     this->throwElementError_m("Marker");
}


inline void ThickTracker::visitMonitor(const Monitor &/*mon*/) {
//     itsOpalBeamline_m.visit(mon, *this, itsBunch_m);
    this->throwElementError_m("Monitor");
}


inline void ThickTracker::visitMultipole(const Multipole &mult) {
    itsOpalBeamline_m.visit(mult, *this, itsBunch_m);

    std::size_t nSlices = mult.getNSlices();
    double length       = mult.getElementLength();
    double gamma        = itsReference.getGamma();
    double p0           = itsReference.getP();
    double q            = itsReference.getQ(); // particle change [e]
    double gradB        = mult.getField().getNormalComponent(2) * ( Physics::c/ p0 ); // [T / m]
    //FIXME remove the next line
    gradB = std::round(gradB*1e6)*1e-6;

    double k1           = gradB * q *Physics::c / p0; // [1 / m^2]

    elements_m.push_back(std::make_tuple(hamiltonian_m.quadrupole(gamma, q, k1),
                                         nSlices,
                                         length));
}

inline void ThickTracker::visitProbe(const Probe &/*probe*/) {
//     itsOpalBeamline_m.visit(prob, *this, itsBunch_m);
    this->throwElementError_m("Probe");
}


inline void ThickTracker::visitRBend(const RBend &/*bend*/) {
//     itsOpalBeamline_m.visit(bend, *this, itsBunch_m);
    this->throwElementError_m("RBend");
}


inline void ThickTracker::visitRFCavity(const RFCavity &/*as*/) {
//     itsOpalBeamline_m.visit(as, *this, itsBunch_m);
    this->throwElementError_m("RFCavity");
}

inline void ThickTracker::visitTravelingWave(const TravelingWave &/*as*/) {
//     itsOpalBeamline_m.visit(as, *this, itsBunch_m);
    this->throwElementError_m("TravelingWave");
}


inline void ThickTracker::visitRFQuadrupole(const RFQuadrupole &/*rfq*/) {
//     itsOpalBeamline_m.visit(rfq, *this, itsBunch_m);
    this->throwElementError_m("RFQuadrupole");
}

inline void ThickTracker::visitSBend(const SBend &bend) {
    itsOpalBeamline_m.visit(bend, *this, itsBunch_m);

    double q     = itsReference.getQ(); // particle change [e]
    double ekin = bend.getDesignEnergy();

    double m     = itsReference.getM(); // eV / c^2
    double gamma = ekin / m + 1.0;
    double beta  = std::sqrt(1.0 - 1.0 / ( gamma * gamma ) );
    double p0    = gamma * beta * m; // eV / c

    double B     = bend.getB() * Physics::c / p0; // T = V * s / m^2
    double r     = std::abs(p0   / ( q * B * Physics::c )); // m

    double k0    = B * q * Physics::c / p0; // V * s * e * m / (m^2 * s * c )

    // [1/m]
    double h = 1.0 / r;

    double L = bend.getElementLength();

    if ( k0 < 0.0 )
        h *= -1.0;

    std::size_t nSlices = bend.getNSlices();

    double arclength    = 2.0 * r * std::asin( L / ( 2.0 * r ) ); //bend.getEffectiveLength();

    // Fringe fields currently not working
    //FIXME e1 not initialised
    //insert Entrance Fringefield
    double e1 = bend.getEntranceAngle();
    if (std::abs(e1) > 1e-6){
        elements_m.push_back(std::make_tuple(hamiltonian_m.fringeField(e1, h),
                                             1, 0.0));
    }

    //insert Dipole "body"
    elements_m.push_back(std::make_tuple(hamiltonian_m.sbend(gamma, h, k0),
                                         nSlices,
                                         arclength));

    //FIXME e2 not initialised
    //insert Exit Fringe field
    double e2 = bend.getExitAngle();
    if (std::abs(e2) > 1e-6){
        elements_m.push_back(std::make_tuple(hamiltonian_m.fringeField(e2, h),
                                             1, 0.0));
    }

}


inline void ThickTracker::visitSeparator(const Separator &/*sep*/) {
//     itsOpalBeamline_m.visit(sep, *this, itsBunch_m);
    this->throwElementError_m("Separator");
}


inline void ThickTracker::visitSeptum(const Septum &/*sept*/) {
//     itsOpalBeamline_m.visit(sept, *this, itsBunch_m);
    this->throwElementError_m("Septum");
}


inline void ThickTracker::visitSolenoid(const Solenoid &/*solenoid*/) {
//     itsOpalBeamline_m.visit(solenoid, *this, itsBunch_m);
    this->throwElementError_m("Solenoid");
}


inline void ThickTracker::visitParallelPlate(const ParallelPlate &/*pplate*/) {
//     itsOpalBeamline_m.visit(pplate, *this, itsBunch_m);
    this->throwElementError_m("ParallelPlate");
}

#endif // OPAL_ThickTracker_HH