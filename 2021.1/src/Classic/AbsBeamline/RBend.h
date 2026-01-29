//
// Class RBend
//   Interface for a rectangular bend magnet.
//
//   A rectangular bend magnet physically has a rectangular shape.
//
//   The standard rectangular magnet, for purposes of definitions, has a field in
//   the y direction. This produces a bend in the horizontal (x) plane. Bends in
//   other planes can be accomplished by rotating the magnet about the axes.
//
//   A positive bend angle is defined as one that bends a beam to the right when
//   looking down (in the negative y direction) so that the beam is bent in the
//   negative x direction. (This definition of a positive bend is the same whether
//   the charge is positive or negative.)
//
//   A zero degree entrance edge angle is parallel to the x direction in an x/y/s
//   coordinate system. A positive entrance edge angle is defined as one that
//   rotates the positive edge (in x) of the angle toward the positive s axis.
//
//   Since the magnet geometry is a fixed rectangle, the exit edge angle is
//   defined by the bend angle of the magnet and the entrance edge angle. In
//   general, the exit edge angle is equal to the bend angle minus the entrance
//   edge angle.
//
//   ------------------------------------------------------------------------
//
//   This class defines two interfaces:
//
//   1) Interface for rectangular magnets for OPAL-MAP.
//
//    Here we specify multipole components about the curved magnet trajectory.
//
//   2) Interface for rectangular magnets for OPAL-T.
//
//   Here we defined the magnet as a field map.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef CLASSIC_RBend_HH
#define CLASSIC_RBend_HH

#include "AbsBeamline/Bend2D.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"

class RBend: public Bend2D {

public:

    /// Constructor with given name.
    explicit RBend(const std::string &name);

    RBend();
    RBend(const RBend &);
    virtual ~RBend();

    /// Apply visitor to RBend.
    virtual void accept(BeamlineVisitor &) const override;

    /*
     * Methods for OPAL-MAP
     * ====================
     */

    /// Get dipole field of RBend.
    virtual double getB() const = 0;

    /// Get RBend geometry.
    //  Version for non-constant object.
    virtual RBendGeometry &getGeometry() override = 0;

    /// Get RBend geometry
    //  Version for constant object.
    virtual const RBendGeometry &getGeometry() const override = 0;

    /// Get multipole expansion of field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField() override = 0;

    /// Get multipole expansion of field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const override = 0;

    /// Get normal component.
    //  Return the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getNormalComponent(int) const;

    /// Get skew component.
    //  Return the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getSkewComponent(int) const;

    /// Set normal component.
    //  Set the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setNormalComponent(int, double);

    /// Set skew component.
    //  Set the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setSkewComponent(int, double);

    /// Get pole entry face rotation.
    //  Return the rotation of the entry pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getEntryFaceRotation() const = 0;

    /// Get exit pole face rotation.
    //  Return the rotation of the exit pole face with respect to the x-axis.
    //  A positive angle rotates the pole face normal away from the centre
    //  of the machine.
    virtual double getExitFaceRotation() const = 0;

    /// Get entry pole face curvature.
    //  Return the curvature of the entry pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getEntryFaceCurvature() const = 0;

    /// Get exit pole face curvature.
    //  Return the curvature of the exit pole face.
    //  A positive curvature creates a convex pole face.
    virtual double getExitFaceCurvature() const = 0;

    /// Get number of slices.
    //  Slices and stepsize used to determine integration step.
    virtual double getSlices() const = 0;

    /// Get stepsize.
    //  Slices and stepsize used to determine integration step.
    virtual double getStepsize() const = 0;

    /*
     * Methods for OPAL-T.
     */

    virtual ElementBase::ElementType getType() const override;
    virtual void setBendAngle(double angle) override;
    virtual void setEntranceAngle(double entranceAngle) override;

private:

    // Not implemented.
    void operator=(const RBend &);

    virtual bool findChordLength(double &chordLength) override;
};

#endif // CLASSIC_RBend_HH