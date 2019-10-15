#ifndef OPAL_ThickMapper_HH
#define OPAL_ThickMapper_HH

// ------------------------------------------------------------------------
// $RCSfile: ThickMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickMapper
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Mapper.h"

class BMultipoleField;
class PlanarArcGeometry;


// Class ThickMapper
// ------------------------------------------------------------------------
/// Build a map using a finite-length lens for each element.
// Multipole-like elements are done by expanding the Lie series.
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
// [li] All elements are represented by maps for finite-length elements.
//   For multipole-like elements the Lie series is used.
// [li] Geometric transformations ignore rotations about transverse axes and
//   translations along the design orbit and truncate after second order.
// [li] Beam-beam elements are two-dimensional, and the second moment <x,y>
//   of the opposite bunches vanish.
// [/ul]
//
// On going through an element, we increment the map using the following steps:
// Remove the closed orbit from the map.
// Construct the Hamiltonian H about that closed orbit.
//   To do this properly, we must expand the square-root AFTER we translate
//   its argument.  In addition, note that the vector potential will be
//   computed to an order determined by the number of field coefficients
//   handed to it, and this may or may not agree with the order of the map.
//   So as to avoid artificially truncating our map, we make the vector
//   potential exact.  Then when we add it to the Hamiltonian, after
//   translating it according to the closed orbit, the Hamiltonian's
//   truncation order will remain unchanged.
// Compute the map exp(-l:H:)z for the present element.
// Propagate the map by sustituting it into the map for the present element.
//   exp(:f:) exp(:g:)z = G(F(z))
//     |||      |||
//   curr.map  element
// To complete the map, we propagate the closed orbit and add that to the map.

class ThickMapper: public Mapper {

public:

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    explicit ThickMapper(const Beamline &bl, const PartData &data,
                         bool backBeam, bool backTrack);

    virtual ~ThickMapper();

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

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);

private:

    // Not implemented.
    ThickMapper();
    ThickMapper(const ThickMapper &);
    void operator=(const ThickMapper &);

    // Apply drift length.
    void applyDrift(double length);

    // Fringe fields for entrance and exit of a magnetic element.
    void applyEntranceFringe(double edge, double curve,
                             const BMultipoleField &field, double scale);
    void applyExitFringe(double edge, double curve,
                         const BMultipoleField &field, double scale);

};

#endif // OPAL_ThickMapper_HH