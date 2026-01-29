#ifndef CLASSIC_StraightGeometry_HH
#define CLASSIC_StraightGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: StraightGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StraightGeometry
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"

// Class StraightGeometry
// ------------------------------------------------------------------------
/// A geometry representing a straight line.
//  StraightGeometry represents a straight line segment along the local
//  z-axis with no torsion. The origin is defined at the centre of the
//  segment, and the geometry runs from -length/2 to +length/2.  All
//  transformations are correspondingly only simple translations along
//  the z-axis.

class StraightGeometry : public BGeometryBase {
public:

    /// Constructor.
    //  Using the [b]length[/b].
    StraightGeometry(double length = 0.0);

    StraightGeometry(const StraightGeometry &right);
    virtual ~StraightGeometry();
    const StraightGeometry &operator=(const StraightGeometry &right);

    /// Get arc length.
    //  Return the arc length, identical to the straight design length.
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the straight design length.
    virtual double getElementLength() const;

    /// Set design length.
    //  Assign the straight design length.
    virtual void setElementLength(double length);

    /// Get origin.
    //  Return the arc length from the entrance to the origin of the element
    //  (origin >= 0)
    double getOrigin() const;

    /// Get entrance.
    //  Return the arc length from the origin to the entrance of the element
    //  (entrance <= 0)
    double getEntrance() const;

    /// Get exit.
    //  Return the arc length from the origin to the exit of the element
    //  (exit >= 0)
    double getExit() const;

    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the local coordinate system from the
    //  entrance to the exit of the element.
    Euclid3D getTotalTransform() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  origin to the entrance of the element.
    Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  origin to the exit of the element.
    Euclid3D getExitFrame()     const;

private:

    // The design length.
    double len;
};

// inlined (trivial) member functions

inline StraightGeometry::StraightGeometry(double l): len(l)
{}

inline StraightGeometry::StraightGeometry(const StraightGeometry &right):
    BGeometryBase(right), len(right.len)
{}

inline const StraightGeometry &StraightGeometry::operator=
(const StraightGeometry &rhs) {
    len = rhs.len;
    return *this;
}

#endif // CLASSIC_StraightGeometry_HH
