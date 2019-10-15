#ifndef CLASSIC_BeamStrippingRep_HH
#define CLASSIC_BeamStrippingRep_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamStrippingRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamStrippingRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamStripping.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/NullField.h"


// Class BeamStrippingRep
// ------------------------------------------------------------------------
/// Representation for a beam stripping.

class BeamStrippingRep: public BeamStripping {

public:

    /// Constructor with given name.
    explicit BeamStrippingRep(const std::string &name);

    BeamStrippingRep();
    BeamStrippingRep(const BeamStrippingRep &);
    virtual ~BeamStrippingRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

//    virtual double getPressure() const;
//    void setPressure(double);

    /// Get field.
    //  Version for non-constant object.
    virtual NullField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const NullField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual PlanarArcGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const PlanarArcGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;


private:

    // Not implemented.
    void operator=(const BeamStrippingRep &);

    // The zero magnetic field.
    NullField field;

    // The geometry.
    PlanarArcGeometry geometry;

};

#endif // CLASSIC_BeamStrippingRep_HH