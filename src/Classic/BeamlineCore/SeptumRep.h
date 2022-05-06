//
// Class SeptumRep
//   Representation for Septum.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef CLASSIC_SeptumRep_HH
#define CLASSIC_SeptumRep_HH

#include "AbsBeamline/Septum.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

class SeptumRep: public Septum {

public:

    /// Constructor with given name.
    explicit SeptumRep(const std::string &name);

    SeptumRep();
    SeptumRep(const SeptumRep &);
    virtual ~SeptumRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
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
    void operator=(const SeptumRep &);



};

#endif // CLASSIC_SeptumRep_HH
