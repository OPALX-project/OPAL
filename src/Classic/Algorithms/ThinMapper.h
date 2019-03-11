#ifndef CLASSIC_ThinMapper_HH
#define CLASSIC_ThinMapper_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ThinMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThinMapper
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/Mapper.h"

class BMultipoleField;


// Class ThinMapper
// ------------------------------------------------------------------------
/// Construct thin lens map.
//  The visitor class for building a FVps<double,6> for a beamline
//  using a thin-lens approximation for all elements.
//  [P]
//  Approximations used:
//  [UL]
//  [LI]All active elements are represented as thin lenses, sandwiched
//    between two drifts, each half of the element length.
//  [LI]Drifts are handled with a second-order approximation.
//  [LI]Geometric transformations ignore rotations about transverse axes
//    and translations along the design orbit and truncate after second order.
//  [/UL]

class ThinMapper: public Mapper {

public:

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTracl[/b] is true, we track against the beam.
    ThinMapper(const Beamline &bl, const PartData &data,
               bool revBeam, bool revTrack);

    virtual ~ThinMapper();


    /// Apply the algorithm to a beam-beam.
    virtual void visitBeamBeam(const BeamBeam &);

    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to a corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to a drift.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a diagnostic.
    virtual void visitDiagnostic(const Diagnostic &);

    /// Apply the algorithm to a drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a drift.
    virtual void visitUndulator(const Undulator &);

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

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &);

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

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &);

protected:

    // Apply a drift length.
    // Approximate method to gain speed.
    void applyDrift(double length);

private:

    // Not implemented.
    ThinMapper();
    ThinMapper(const ThinMapper &);
    void operator=(const ThinMapper &);
};

#endif // CLASSIC_ThinMapper_HH
