#ifndef CLASSIC_NullGeometry_HH
#define CLASSIC_NullGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: NullGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NullGeometry
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

class Euclid3D;


// Class NullGeometry
// ------------------------------------------------------------------------
/// Geometry representing an identity transform.

class NullGeometry: public BGeometryBase {
public:

    NullGeometry();
    NullGeometry(const NullGeometry &);
    virtual ~NullGeometry();
    const NullGeometry &operator=(const NullGeometry &);

    /// Get arc length.
    //  Always zero for this class.
    virtual double getArcLength() const;

    /// Get design length.
    //  Always zero for this class.
    virtual double getElementLength() const;

    /// Get origin position.
    //  Always zero for this class.
    virtual double getOrigin() const;

    /// Get transform.
    //  transform from [b]fromS[/b] to [b]toS[/b] is always an identity
    //  for this class.
    virtual Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    virtual Euclid3D getTransform(double s) const;
};


// inlined (trivial) member functions
inline NullGeometry::NullGeometry()
{}

inline NullGeometry::NullGeometry(const NullGeometry & o) : BGeometryBase(o)
{}

inline const NullGeometry &
NullGeometry::operator=(const NullGeometry &) {
    return *this;
}

#endif // CLASSIC_NullGeometry_HH
