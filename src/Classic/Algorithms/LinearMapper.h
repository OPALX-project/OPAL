#ifndef MAD_LinearMapper_HH
#define MAD_LinearMapper_HH

// ------------------------------------------------------------------------
// $RCSfile: LinearMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LinearMapper
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/31 07:12:38 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractMapper.h"
#include "FixedAlgebra/LinearMap.h"


class BMultipoleField;
class Euclid3D;

template <class T, int N> class FTps;
template <class T, int N> class FVps;


// Class LinearMapper
// ------------------------------------------------------------------------
/// Build a map using a linear map for each element.
// Phase space coordinates numbering:
// [tab 3 b]
// [row]number [&]name        [&]unit  [/row]
// [row]0      [&]x           [&]metres [/row]
// [row]1      [&]p_x/p_r     [&]1      [/row]
// [row]2      [&]y           [&]metres [/row]
// [row]3      [&]p_y/p_r     [&]1      [/row]
// [row]4      [&]v*delta_t   [&]metres [/row]
// [row]5      [&]delta_p/p_r [&]1      [/row]
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
// [li] All elements are represented by maps for finite-length elements.
// [li] Geometric transformations ignore rotations about transverse axes and
//   translations along the design orbit and truncate after second order.
// [li] Beam-beam elements are two-dimensional, and the second moment <x,y>
//   of the opposite bunches vanish.
// [/ul]

class LinearMapper: public AbstractMapper {

public:

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    LinearMapper(const Beamline &beamline, const PartData &,
                 bool revBeam, bool revTrack);

    virtual ~LinearMapper();

    /// Return the linear part of the accumulated map.
    virtual void getMap(LinearMap<double, 6> &) const;

    /// Return the full map accumulated so far.
    virtual void getMap(FVps<double, 6> &) const;

    /// Reset the linear part of the accumulated map for restart.
    virtual void setMap(const LinearMap<double, 6> &);

    /// Reset the full map for restart.
    virtual void setMap(const FVps<double, 6> &);

    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a beam stripping.
    virtual void visitBeamStripping(const BeamStripping &);


    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to an arbitrary component.
    //  This override calls the component to track the map.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a Corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a drift.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a Diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a Drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a undulator.
    virtual void visitUndulator(/* const Undulator & */);

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

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend0(const RBend &);
    virtual void visitRBend(const RBend &);   // exp(:-H:) version
    /// Apply the algorithm to a RFCavity.
    virtual void visitRFCavity(const RFCavity &);

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

    /// Apply the algorithm to an offset beamline object wrapper.
    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to an integrator capable of mapping.
    virtual void visitMapIntegrator(const MapIntegrator &);

protected:

    /// Apply drift length.
    // Propagate current map through a drift.
    void applyDrift(double length);

    /// Transforms fringing fields.
    void applyEntranceFringe(double edge,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge,
                         const BMultipoleField &field, double scale);

    /// Apply linear map, defined by the linear expansions Fx and Fy.
    void applyLinearMap(double length, double refLength, double h,
                        const FTps<double, 2> &Fx, const FTps<double, 2> &Fy);

    /// Apply body of SBend.
    void applyMultipoleBody(double length, double refLength,
                            const BMultipoleField &field, double scale);

    /// Apply thin multipole kick (integrated over length) to all particles.
    void applySBendBody(double length, double refLength, double h,
                        const BMultipoleField &field, double scale);

    /// Thin multipole kick.
    //  Apply a thin multipole kick (integrated over length) to current map.
    void applyThinMultipole(const BMultipoleField &field, double factor);

    /// Apply transform.
    //  Propagate current map through a geometric transformation.
    void applyTransform(const Euclid3D &, double refLength);

    /// Construct the vector potential for a SBend.
    FTps<double, 2>
    buildSBendVectorPotential(const BMultipoleField &, double h);


    // The linear map being accumulated.
    LinearMap <double, 6> itsMap;

private:

    // Not implemented.
    LinearMapper();
    LinearMapper(const LinearMapper &);
    void operator=(const LinearMapper &);

    /// Helper function for finding first-order coefficients.
    static void makeFocus
    (double k, double L, double &c, double &s, double &d, double &f);

};

#endif // MAD_LinearMapper_HH
