//
// Class BeamStrippingRep
//   Defines a concrete representation for beam stripping.
//
// Copyright (c) 2018-2019, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#ifndef CLASSIC_BeamStrippingRep_HH
#define CLASSIC_BeamStrippingRep_HH

#include "AbsBeamline/BeamStripping.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/NullField.h"

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