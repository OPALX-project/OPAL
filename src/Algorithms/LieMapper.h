#ifndef OPAL_LieMapper_HH
#define OPAL_LieMapper_HH

// ------------------------------------------------------------------------
// $RCSfile: LieMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LieMapper
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Mapper.h"
#include "FixedAlgebra/DragtFinnMap.h"

class BMultipoleField;
class PlanarArcGeometry;


// Class LieMapper
// ------------------------------------------------------------------------
/// Build a Lie-algebraic map using a finite-length lens for each elements.
//  All maps are factored in ascending Dragt-Finn order.
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
// All elements are represented by maps for finite-length elements.
// Several important pieces are still MISSING from this algorithm.

class LieMapper: public AbstractMapper {

public:

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    LieMapper(const Beamline &bl, const PartData &data,
              bool backBeam, bool backTrack, int order);

    virtual ~LieMapper();

    /// Return the linear part of the accumulated map.
    virtual void getMap(LinearMap<double, 6> &) const;

    /// Return the accumulated map.
    virtual void getMap(FVps<double, 6> &) const;

    /// Return the full map accumulated so far.
    virtual void getMap(DragtFinnMap<3> &) const;

    /// Reset the linear part of the accumulated map for restart.
    virtual void setMap(const LinearMap<double, 6> &);

    /// Reset the accumulated map for restart.
    virtual void setMap(const FVps<double, 6> &);

    /// Reset the full map for restart.
    virtual void setMap(const DragtFinnMap<3> &);


    /// Apply the algorithm to a BeamBeam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a BeamStripping.
    virtual void visitBeamStripping(const BeamStripping &);

    /// Apply the algorithm to a Collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to a Component.
    virtual void visitComponent(const Component &);

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

    /// Apply the algorithm to a Patch.
    virtual void visitPatch(const Patch &);

    /// Apply the algorithm to a Probe.
    virtual void visitProbe(const Probe &);

    /// Apply the algorithm to a RBend.
    virtual void visitRBend(const RBend &);

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

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

private:

    // Not implemented.
    LieMapper();
    LieMapper(const LieMapper &);
    void operator=(const LieMapper &);

    // Apply drift length.
    // Propagate current map through a drift.
    void applyDrift(double length);

    // Transforms fringing fields.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

    /// Apply transform.
    //  Propagate current map through a geometric transformation.
    //  Linear approximation for the time being.
    virtual void applyTransform(const Euclid3D &, double refLength = 0.0);

    // The Lie map being accumulated.
    DragtFinnMap<3> itsMap;

    // The desired map order.
    int itsOrder;
};

#endif // OPAL_LieMapper_HH