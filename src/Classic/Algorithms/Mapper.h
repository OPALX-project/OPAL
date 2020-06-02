#ifndef CLASSIC_Mapper_HH
#define CLASSIC_Mapper_HH

// ------------------------------------------------------------------------
// $RCSfile: Mapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Mapper
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractMapper.h"
#include "AbsBeamline/Patch.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "Algorithms/PartData.h"

class AlignWrapper;
class BMultipoleField;
class Euclid3D;
template <class T, int N> class LinearMap;
template <class T, int N> class TransportMap;


// Class Mapper
// ------------------------------------------------------------------------
/// Build transfer map.
//  An abstract visitor class implementing the default behaviour for all
//  visitors capable of tracking a transfer map through a beam line.
//  It implements some default behaviour for such visitors.
//  [P]
//  Phase space coordinates (in this order):
//  [DL]
//  [DT]x:[DD]
//    horizontal displacement (metres).
//  [DT]p_x/p_r:[DD]
//    horizontal canonical momentum (no dimension).
//  [DT]y:[DD]
//    vertical displacement (metres).
//  [DT]p_y/p_r:[DD]
//    vertical canonical momentum (no dimension).
//  [DT]delta_p/p_r:[DD]
//    relative momentum error (no dimension).
//  [DT]v*delta_t:[DD]
//    time difference delta_t w.r.t. the reference frame which moves with
//    uniform velocity
//    [P]
//    v_r = c*beta_r = p_r/m
//    [P]
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

class Mapper: public AbstractMapper {

public:

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    Mapper(const Beamline &bl, const PartData &data,
           bool revBeam, bool revTrack);

    virtual ~Mapper();

    /// Return the linear part of the accumulated map.
    virtual void getMap(LinearMap<double, 6> &) const;

    /// Return the second-order part of the accumulated map.
    virtual void getMap(TransportMap<double, 6> &) const;

    /// Return the full map accumulated so far.
    virtual void getMap(FVps<double, 6> &) const;

    /// Reset the linear part of the accumulated map for restart.
    virtual void setMap(const LinearMap<double, 6> &);

    /// Reset the second-order part of the accumulated map for restart.
    virtual void setMap(const TransportMap<double, 6> &);

    /// Reset the full map for restart.
    virtual void setMap(const FVps<double, 6> &);


    /// Apply the algorithm to an arbitrary component.
    //  This override calls the component to track the map.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat);


    /// Apply the algorithm to an align wrapper.
    virtual void visitAlignWrapper(const AlignWrapper &);

protected:

    /// Apply drift length.
    //  Propagate current map through a drift.
    void applyDrift(double length);

    /// Thin multipole kick.
    //  Apply a thin multipole kick (integrated over length) to current map.
    void applyThinMultipole(const BMultipoleField &field, double factor);

    /// Thin SBend kick.
    //  Special kick routine for thin SBend.
    void applyThinSBend(const BMultipoleField &field, double scale, double h);

    /// Apply transform.
    //  Propagate current map through a geometric transformation.
    void applyTransform(const Euclid3D &, double refLength = 0.0);

    /// The transfer map being built.
    FVps<double, 6> itsMap;

private:

    // Not implemented.
    Mapper();
    Mapper(const Mapper &);
    void operator=(const Mapper &);
};

#endif // CLASSIC_Mapper_HH
