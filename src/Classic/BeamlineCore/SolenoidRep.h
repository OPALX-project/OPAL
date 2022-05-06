//
// Class SolenoidRep
//   Representation for a solenoid magnet.
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
#ifndef CLASSIC_SolenoidRep_HH
#define CLASSIC_SolenoidRep_HH

#include "AbsBeamline/Solenoid.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/ConstBzField.h"


class SolenoidRep: public Solenoid {

public:

    /// Constructor with given name.
    explicit SolenoidRep(const std::string &name);

    SolenoidRep();
    SolenoidRep(const SolenoidRep &);
    virtual ~SolenoidRep();

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
    virtual ConstBzField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const ConstBzField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get field.
    //  Return the solenoid field in Teslas.
    virtual double getBz() const;

    /// Set field.
    //  Assign the solenoid field in Teslas.
    virtual void setBz(double Bz);

private:

    // Not implemented.
    void operator=(const SolenoidRep &);

    /// The solenoid geometry.
    StraightGeometry geometry;

    /// The solenoid field.
    ConstBzField field;
};

#endif // CLASSIC_SolenoidRep_HH
