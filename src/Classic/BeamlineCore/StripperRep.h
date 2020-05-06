//
// Class StripperRep
//   Representation for Stripper
//
// Copyright (c) 2011, Jianjun Yang,
//                     Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#ifndef CLASSIC_StripperRep_HH
#define CLASSIC_StripperRep_HH

#include "AbsBeamline/Stripper.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

class StripperRep: public Stripper {

public:

    /// Constructor with given name.
    explicit StripperRep(const std::string &name);

    StripperRep();
    StripperRep(const StripperRep &);
    virtual ~StripperRep();

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
    //  Return the element geometry.
    //  Version for non-constant object.
    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Set active flag.
    virtual void setActive(bool = true);

protected:

    /// The zero magnetic field.
    NullField field;

    /// The septa's geometry.
    StraightGeometry geometry;

    /// The active/inactive flag.
    bool active;

private:

    // Not implemented.
    void operator=(const StripperRep &);



};

#endif // CLASSIC_StripperRep_HH
