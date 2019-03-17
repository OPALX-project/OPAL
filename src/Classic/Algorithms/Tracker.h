#ifndef CLASSIC_Tracker_HH
#define CLASSIC_Tracker_HH

// ------------------------------------------------------------------------
// $RCSfile: Tracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Tracker
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractTracker.h"
#include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartData.h"
#include "FixedAlgebra/FTps.h"

#include "Utilities/ClassicField.h"

class BMultipoleField;
class Euclid3D;
class OpalParticle;


// Class Tracker
// ------------------------------------------------------------------------
/// Track particles or bunches.
//  An abstract base class for all visitors capable of tracking particles
//  through a beam element.
//  [P]
//  Phase space coordinates (in this order):
//  [DL]
//  [DT]x:[DD]
//    horizontal displacement (metres).
//  [DT]p_x/p_r:[DD]
//     horizontal canonical momentum (no dimension).
//  [DT]y:[DD]
//    vertical displacement (metres).
//  [DT]p_y/p_r:[DD]
//    vertical canonical momentum (no dimension).
//  [DT]delta_p/p_r:[DD]
//    relative momentum error (no dimension).
//  [DT]v*delta_t:[DD]
//    time difference delta_t w.r.t. the reference frame which moves with
//    uniform velocity
//  [P]
//    v_r = c*beta_r = p_r/m
//  [P]
//    along the design orbit, multiplied by the instantaneous velocity v of
//    the particle (metres).
//  [/DL]
//  Where
//  [DL]
//  [DT]p_r:[DD]
//    is the constant reference momentum defining the reference frame velocity.
//  [DT]m:[DD]
//    is the rest mass of the particles.
//  [/DL]
//  Other units used:
//  [DL]
//  [DT]reference momentum:[DD]
//    electron-volts.
//  [DT]accelerating voltage:[DD]
//    volts.
//  [DT]separator voltage:[DD]
//    volts.
//  [DT]frequencies:[DD]
//    hertz.
//  [DT]phase lags:[DD]
//    multiples of (2*pi).
//  [/DL]


class Tracker: public AbstractTracker {

public:

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  The particle bunch is initially empty.
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    Tracker(const Beamline &, const PartData &,
            bool backBeam, bool backTrack);

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  The particle bunch is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    Tracker(const Beamline &, PartBunchBase<double, 3> *bunch,
            const PartData &, bool backBeam, bool backTrack);

    virtual ~Tracker();

    /// Return the current bunch.
    const PartBunchBase<double, 3> *getBunch() const;

    /// Add particle to bunch.
    void addToBunch(const OpalParticle &);

    /// Store the bunch.
    //~ void setBunch(const PartBunch &);

    /// Apply the algorithm to an arbitrary component.
    //  This override calls the component to track the bunch.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat);

    /// Apply the algorithm to an align wrapper.
    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply the algorithm to an integrator capable of tracking.
    virtual void visitTrackIntegrator(const TrackIntegrator &);

    /// Apply the algorithm to an integrator capable of mapping.
    virtual void visitMapIntegrator(const MapIntegrator &);

    /// set total number of tracked bunches
    virtual void setNumBunch(int) {};

    /// get total number of tracked bunches
    virtual int  getNumBunch() { return 0; }

    /// set multipacting flag
    //virtual void setMpacflg(bool /* mpacflg */) {};

    // standing wave structures
    FieldList cavities_m;

    const Beamline &itsBeamline_m;

protected:

    /// Apply a drift length.
    void applyDrift(double length);

    // Apply thin multipole kick.                                         .
    void applyThinMultipole(const BMultipoleField &field, double factor);

    // Special kick routine for thin SBend.
    void applyThinSBend(const BMultipoleField &field, double scale, double h);

    /// Apply a geometric transformation.
    void applyTransform(const Euclid3D &, double refLength = 0.0);

    /// Construct vector potential for a Multipole.
    FTps<double, 2> buildMultipoleVectorPotential2D(const BMultipoleField &);

    /// Construct vector potential for a Multipole.
    FTps<double, 6> buildMultipoleVectorPotential(const BMultipoleField &);

    /// Construct vector potential for a SBend.
    FTps<double, 2> buildSBendVectorPotential2D(const BMultipoleField &, double h);

    /// Construct vector potential for a SBend.
    FTps<double, 6> buildSBendVectorPotential(const BMultipoleField &, double h);

    /// The bunch of particles to be tracked.
    PartBunchBase<double, 3>* itsBunch_m;
    //  typedef PartBunch::iterator iterator;

private:

    // Not implemented.
    Tracker();
    Tracker(const Tracker &);
    void operator=(const Tracker &);
};

#endif // CLASSIC_Tracker_HH
