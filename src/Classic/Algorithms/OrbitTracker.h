#ifndef CLASSIC_OrbitTracker_HH
#define CLASSIC_OrbitTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: OrbitTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OrbitTracker
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractTracker.h"
#include "FixedAlgebra/FVector.h"

class BMultipoleField;
class Euclid3D;

template <class T, int N> class FTps;


// Class OrbitTracker
// ------------------------------------------------------------------------
/// Track closed orbit.
//  Uses a thick lens method to find the orbit for each element.

class OrbitTracker: public AbstractTracker {

public:

    // Particle coordinate numbers.
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    OrbitTracker(const Beamline &, const PartData &,
                 bool backBeam, bool backTrack);

    virtual ~OrbitTracker();


    /// Return the current orbit.
    const FVector<double, 6> &getOrbit() const;

    /// Reset the current orbit.
    void setOrbit(const FVector<double, 6> orbit);

    /// Apply the algorithm to a beam-beam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a beam stripping.
    virtual void visitBeamStripping(const BeamStripping &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a drift.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a undulator - NOT implemented.
    virtual void visitUndulator(/* const Undulator & */);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a Lambertson.
    virtual void visitLambertson(const Lambertson &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a RF quadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a separator.
    virtual void visitSeparator(const Separator &);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &);

    /// Apply the algorithm to an align wrapper..
    virtual void visitAlignWrapper(const AlignWrapper &);

private:

    // Not implemented.
    OrbitTracker();
    OrbitTracker(const OrbitTracker &);
    void operator=(const OrbitTracker &);

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

    /// Helper function for finding first-order coefficients.
    static void makeFocus
    (double k, double L, double &c, double &s, double &d, double &f);


    // The current orbit.
    FVector<double, 6> itsOrbit;
};

#endif // CLASSIC_OrbitTracker_HH
