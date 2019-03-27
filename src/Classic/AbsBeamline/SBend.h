#ifndef CLASSIC_SBend_HH
#define CLASSIC_SBend_HH

// ------------------------------------------------------------------------
// $RCSfile: SBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Bend.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"
#include <string>

/*
 * Class SBend
 *
 * Interface for sector bend magnet.
 *
 * A sector bend magnet has a curved geometry. A sector magnet with zero degree
 * edge angles is simply a section of a circle when projected onto the y axis.
 *
 * The standard sector magnet, for purposes of definitions, has a field in the y
 * direction. This produces a bend in the horizontal (x) plane. Bends in other
 * planes can be accomplished by rotating the magnet about the axes.
 *
 * A positive bend angle is defined as one that bends a beam to the right when
 * looking down (in the negative y direction) so that the beam is bent in the
 * negative x direction. (This definition of a positive bend is the same whether
 * the charge is positive or negative.)
 *
 * A zero degree entrance edge angle is parallel to the x direction in an x/y/s
 * coordinate system. A positive entrance edge angle is defined as one that
 * rotates the positive edge (in x) of the angle toward the positive s axis.
 *
 * A zero degree exit edge angle is parallel to the x direction in an x/y/s
 * coordinate system. A positive exit edge angle is defined as one that rotates
 * the positive edge (in x) of the angle toward the negative s axis.
 *
 * ------------------------------------------------------------------------
 *
 * This class defines two interfaces:
 *
 * 1) Interface for sector magnets for OPAL-MAP.
 *
 *  Here we specify multipole components about the curved magnet trajectory.
 *
 *
 * 2) Interface for sector magnets for OPAL-T.
 *
 * Here we defined the magnet as a field map.
 */

class SBend: public Bend {

public:

    /// Constructor with given name.
    explicit SBend(const std::string &name);

    SBend();
    SBend(const SBend &);
    virtual ~SBend();

    /// Apply visitor to SBend.
    virtual void accept(BeamlineVisitor &) const override;


    /*
     * Methods for OPAL-MAP
     * ====================
     */

    /// Get dipole field of SBend.
    virtual double getB() const = 0;

    /// Get SBend geometry.
    //  Version for non-constant object.
    virtual PlanarArcGeometry &getGeometry() override = 0;

    /// Get SBend geometry
    //  Version for constant object.
    virtual const PlanarArcGeometry &getGeometry() const override = 0;

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
     * Methods for OPAL-SLICE.
     */
    virtual void addKR(int /* i */, double /* t */, Vector_t &/* K */) override { };
    virtual void addKT(int /* i */, double /* t */, Vector_t &/* K */) override { };


    virtual ElementBase::ElementType getType() const override;

private:

    // Not implemented.
    void operator=(const SBend &);

    virtual bool findChordLength(Inform &msg,
                                 double &chordLength) override;

};

#endif // CLASSIC_SBend_HH
