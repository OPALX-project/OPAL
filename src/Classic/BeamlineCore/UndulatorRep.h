#ifndef CLASSIC_UndulatorRep_HH
#define CLASSIC_UndulatorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: UndulatorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UndulatorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Undulator.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"


// Class UndulatorRep
// ------------------------------------------------------------------------
/// Representation for a drift space.

class UndulatorRep: public Undulator {

public:

    /// Constructor with given name.
    explicit UndulatorRep(const std::string &name);

    UndulatorRep();
    UndulatorRep(const UndulatorRep &);
    virtual ~UndulatorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual NullField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const NullField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

private:

    // Not implemented.
    void operator=(const UndulatorRep &);

    /// The zero magnetic field.
    NullField field;

    /// The geometry.
    StraightGeometry geometry;
};

#endif // CLASSIC_UndulatorRep_HH
