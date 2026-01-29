#ifndef CLASSIC_BeamlineGeometry_HH
#define CLASSIC_BeamlineGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamlineGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamlineGeometry
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"

class Beamline;


// Class BeamlineGeometry
// ------------------------------------------------------------------------
/// Implements the composite geometry of a beam line.

class BeamlineGeometry: public BGeometryBase {

public:

    /// Constructor.
    //  The geometry is linked to th beamline [b]line[/b].
    explicit BeamlineGeometry(const Beamline &line);

    virtual ~BeamlineGeometry();

    /// Get arc length.
    //  Return the length of the geometry, measured along the design orbit.
    virtual double getArcLength() const;

    /// Get element length.
    //  Return the length of the geometry, measured along the design polygone.
    virtual double getElementLength() const;

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the
    //  geometry (non-negative).
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the
    //  geometry (non-positive0)
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the
    //  geometry (non-negative).
    virtual double getExit() const;

    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    virtual Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    virtual Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the local coordinate system from the
    //  entrance to the exit of the element.
    virtual Euclid3D getTotalTransform() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  origin to the entrance of the element.
    virtual Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  origin to the exit of the element.
    virtual Euclid3D getExitFrame() const;

private:

    // Not implemented.
    BeamlineGeometry();
    BeamlineGeometry(const BeamlineGeometry &);
    void operator=(const BeamlineGeometry &);

    // The beamline whose geometry [b]this[/b] represents.
    const Beamline &itsLine;
};

#endif // CLASSIC_BeamlineGeometry_HH
