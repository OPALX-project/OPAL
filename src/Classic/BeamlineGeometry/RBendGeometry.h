#ifndef CLASSIC_RBendGeometry_HH
#define CLASSIC_RBendGeometry_HH

// ------------------------------------------------------------------------
// $RCSfile: RBendGeometry.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendGeometry
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/StraightGeometry.h"

// Class RBendGeometry
// ------------------------------------------------------------------------
/// The geometry for a RBend element.
//  It is a Geometry wrapper which adds two rotations by alpha/2 (local
//  y-axis rotations) to the entrance and exit planes of a StraightGeometry.
//  The y-rotations become part of the global geometry definition.
//  {P}
//  NOTE: in general the transformations returned include the effects of
//  the y-rotations when the required point specifies the entrance or exit
//  point. Requests for transformations within the geometry (i.e. from s1
//  to s2, where s1 and/or s2 are not the entrance or exit planes) do not
//  contain the y-rotations.
//  {P}
//  A RBendGeometry can be seen as an OffsetGeometry, whose global geometry
//  is a PlanarArcGeometry, and whose local geometry is a StraightGeometry.

class RBendGeometry : public StraightGeometry {
public:

    /// Constructor.
    //  Construct an RBendGeometry from [b]length[/b] and [b]angle[/b].
    RBendGeometry(double length, double angle);

    RBendGeometry(const RBendGeometry &);
    virtual ~RBendGeometry();
    const RBendGeometry &operator=(const RBendGeometry &);

    /// Get arc length.
    //  Return the length measured along a circular arc tangent to the
    //  local s-axis at entrance and exit; an approximation to the actual
    //  design orbit.
    virtual double getArcLength() const;

    /// Get angle.
    //  Return the total bend angle.
    virtual double getBendAngle() const;

    /// Set angle.
    //  Assign the bend angle.
    void setBendAngle(double angle);

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
    Euclid3D getExitFrame() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the global geometry to the local geometry at entrance.
    Euclid3D getEntrancePatch() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the local geometry to the global geometry at exit.
    Euclid3D getExitPatch()     const;

private:

    double half_angle;
};

#endif // CLASSIC_RBendGeometry_HH