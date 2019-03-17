#ifndef CLASSIC_CyclotronWrapper_HH
#define CLASSIC_CyclotronWrapper_HH

// ------------------------------------------------------------------------
// $RCSfile: CyclotronWrapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronWrapper
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Cyclotron.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"
#include "MemoryManagement/Pointer.h"

class BMultipoleField;

// Class CyclotronWrapper
// ------------------------------------------------------------------------
/// Representation of a perturbed cyclotron.
//  A CyclotronWrapper represents a unique instance of a cyclotron
//  in the accelerator model. It defines imperfections of the field,
//  related to an ``ideal'' magnet contained in the wrapper.

class CyclotronWrapper: public Cyclotron {

public:

    /// Constructor.
    //  Constructs a wrapper for its argument.
    //  The wrapper is not sharable by default.
    explicit CyclotronWrapper(Cyclotron *);

    CyclotronWrapper(const CyclotronWrapper &);
    virtual ~CyclotronWrapper();


    /// Apply visitor to modified bend.
    virtual void accept(BeamlineVisitor &) const;

    /// Make clone.
    virtual ElementBase *clone() const;

    /// Make structural copy.
    virtual ElementBase *copyStructure();

    /// Get multipole field error.
    //  This method can be used to get or set the error field. The error field
    //  offset is mutable, so as to allow changing it in a constant structure.
    virtual BMultipoleField &errorField() const ;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Get field.
    //  Return the perturbed field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Return the perturbed field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual PlanarArcGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const PlanarArcGeometry &getGeometry() const;

    /// Get number of slices.
    virtual double getSlices() const;

    /// Get stepsize.
    virtual double getStepsize() const;

    /// Get element type string.
    virtual ElementBase::ElementType getType() const;

    /// Get design Cyclotron.
    //  Version for constant object.
    virtual const Cyclotron &getDesign() const;

    /// Get design Cyclotron.
    //  Version for non-constant object.
    virtual Cyclotron &getDesign();

    /// Make wrapper for this bend.
    //  Return [b]this[/b], since this is already a field wrapper.
    virtual ElementBase *makeFieldWrapper();

    /// Remove field wrapper.
    virtual ElementBase *removeFieldWrapper();

    /// Remove field wrapper.
    virtual const ElementBase *removeFieldWrapper() const;

    /// Remove all wrappers.
    virtual ElementBase *removeWrappers();

    /// Remove all wrappers.
    virtual const ElementBase *removeWrappers() const;

private:

    // Not implemented.
    CyclotronWrapper();
    void operator=(const CyclotronWrapper &);

    // The ``design'' magnet.
    Pointer<Cyclotron> itsDesign;

    /// The error field.
    mutable BMultipoleField itsError;

    // The magnetic field after application of all modifiers.
    mutable BMultipoleField tempField;
};

#endif // CLASSIC_CyclotronWrapper_HH