#ifndef CLASSIC_AlignWrapper_HH
#define CLASSIC_AlignWrapper_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignWrapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignWrapper
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "MemoryManagement/Pointer.h"

class BeamlineVisitor;


// Class AlignWrapper
// ------------------------------------------------------------------------
/// Define the position of a misaligned element.
//  An AlignWrapper is used to store misalignment errors or deliberate
//  misalignments.  It acts as a wrapper around a component or a complete
//  beam line.  Rotations and translations are defined about the design
//  local frame, which in turn is specified by the position of the element
//  on the design geometry.  An AlignWrapper is non-sharable by default.

class AlignWrapper: public ElementBase {

    friend class ElementBase;

public:

    /// Apply BeamlineVisitor.
    virtual void accept(BeamlineVisitor &) const override;

    /// Return clone.
    //  Return an identical deep copy of the wrapper and its contents.
    virtual AlignWrapper *clone() const override;

    /// Make structural copy.
    virtual ElementBase *copyStructure() override;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable() override;


    /// Get entrance patch.
    //  Returns the entrance patch (transformation) which is used to
    //  transform the global geometry to the local geometry at entrance
    //  of the misaligned element.
    virtual Euclid3D getEntranceTransform() const;

    /// Get exit patch.
    //  Returns the exit patch (transformation) which is used to
    //  transform the local geometry to the global geometry at exit of
    //  the misaligned element.
    virtual Euclid3D getExitTransform() const;


    /// Return the contained element.
    virtual ElementBase *getElement() const;

    /// Replace the contained element.
    void setElement(ElementBase *);

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual BGeometryBase &getGeometry() override;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for constant object.
    virtual const BGeometryBase &getGeometry() const override;

    /// Return the offset.
    //  This method can be used to get or set the offset. The offset is
    //  declared as mutable, so as to allow changing it in a constant
    //  structure.
    Euclid3D &offset() const;

    /// Get element type std::string.
    //  Returns the type std::string for the enclosed item.
    virtual ElementBase::ElementBase::ElementType getType() const override;

    /// Disallow misalignment of an already misaligned object.
    //  This method returns [b]this[/b], since "this" is already an
    //  AlignWrapper.
    virtual ElementBase *makeAlignWrapper() override;

    /// Allow field errors.
    //  Wrap the contained element in a field wrapper, unless such a
    //  wrapper already exists.
    virtual ElementBase *makeFieldWrapper() override;

    /// Remove AlignWrapper.
    //  Return the element or field wrapper contained in "this".
    //  Version for non-const object.
    virtual ElementBase *removeAlignWrapper() override;

    /// Remove AlignWrapper.
    //  Return the element or field wrapper contained in "this".
    //  Version for const object.
    virtual const ElementBase *removeAlignWrapper() const override;

    /// Remove field wrapper.
    //  Remove any field wrapper on the contained object.
    virtual ElementBase *removeFieldWrapper() override;

    /// Remove field wrapper.
    //  Remove the field wrapper for constant object.
    virtual const ElementBase *removeFieldWrapper() const override;

    /// Return the design element.
    //  This method removes all wrappers on the contained object.
    //  Version for non-const object.
    virtual ElementBase *removeWrappers() override;

    /// Return the design element.
    //  Version for const object.
    virtual const ElementBase *removeWrappers() const override;

private:

    AlignWrapper(const AlignWrapper &);
    virtual ~AlignWrapper();

    // Constructor.
    // Constructs an AlignWrapper wrapping the given element.
    AlignWrapper(ElementBase *);

    // Not implemented.
    AlignWrapper();
    void operator=(const AlignWrapper &);

    // The pointer to the wrapped element.
    Pointer<ElementBase> itsElement;

    // The offset of the contained element.
    // The data member [b]offset[/b] is declared mutable, so as to allow
    // setting of a misalignment by a Visitor (algorithm) in an otherwise
    // constant line.  This mechanism protects the beam line structure against
    // accidental change by a Visitor.
    mutable Euclid3D itsOffset;
};

#endif // CLASSIC_AlignWrapper_HH
